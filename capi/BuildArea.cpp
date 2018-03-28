/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 * 
 * Ported from rtgeom_geos.c from
 *   rttopo - topology library
 *   http://git.osgeo.org/gitea/rttopo/librttopo
 * with relicensing from GPL to LGPL with Copyright holder permission.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Copyright 2011-2014 Sandro Santilli <strk@kbt.io>
 *
 **********************************************************************/

#include "geos_c.h"

#include <geos/util/GEOSContextHandle.h>

/* ------------ BuildArea stuff ---------------------------------------------------------------------{ */

typedef struct Face_t {
  const GEOSGeometry* geom;
  GEOSGeometry* env;
  double envarea;
  struct Face_t* parent; /* if this face is an hole of another one, or NULL */
} Face;

static Face* newFace(GEOSContextHandle_t handle, const GEOSGeometry* g);
static void delFace(GEOSContextHandle_t handle, Face* f);
static unsigned int countParens(GEOSContextHandle_t handle, const Face* f);
static void findFaceHoles(GEOSContextHandle_t handle, Face** faces, int nfaces);

static Face*
newFace(GEOSContextHandle_t handle, const GEOSGeometry* g)
{
  Face* f = new Face;
  f->geom = g;
  f->env = GEOSEnvelope_r(handle, f->geom);
  GEOSArea_r(handle, f->env, &f->envarea);
  f->parent = NULL;
  /* rtnotice(ctx, "Built Face with area %g and %d holes", f->envarea, GEOSGetNumInteriorRings_r(handle, f->geom)); */
  return f;
}

static unsigned int
countParens(GEOSContextHandle_t, const Face* f)
{
  unsigned int pcount = 0;
  while ( f->parent ) {
    ++pcount;
    f = f->parent;
  }
  return pcount;
}

/* Destroy the face and release memory associated with it */
static void
delFace(GEOSContextHandle_t handle, Face* f)
{
  GEOSGeom_destroy_r(handle, f->env);
  delete f;
}


static int
compare_by_envarea(const void* g1, const void* g2)
{
  Face* f1 = *(Face**)g1;
  Face* f2 = *(Face**)g2;
  double n1 = f1->envarea;
  double n2 = f2->envarea;

  if ( n1 < n2 ) return 1;
  if ( n1 > n2 ) return -1;
  return 0;
}

/* Find holes of each face */
static void
findFaceHoles(GEOSContextHandle_t handle, Face** faces, int nfaces)
{
  int i, j, h;

  /* We sort by envelope area so that we know holes are only
   * after their shells */
  qsort(faces, nfaces, sizeof(Face*), compare_by_envarea);
  for (i=0; i<nfaces; ++i) {
    Face* f = faces[i];
    int nholes = GEOSGetNumInteriorRings_r(handle, f->geom);
#ifdef RTDEBUGF
    RTDEBUGF(ctx, 2, "Scanning face %d with env area %g and %d holes", i, f->envarea, nholes);
#endif
    for (h=0; h<nholes; ++h) {
      const GEOSGeometry *hole = GEOSGetInteriorRingN_r(handle, f->geom, h);
#ifdef RTDEBUGF
      RTDEBUGF(ctx, 2, "Looking for hole %d/%d of face %d among %d other faces", h+1, nholes, i, nfaces-i-1);
#endif
      for (j=i+1; j<nfaces; ++j) {
        const GEOSGeometry *f2er;
        Face* f2 = faces[j];
        if ( f2->parent ) continue; /* hole already assigned */
        f2er = GEOSGetExteriorRing_r(handle, f2->geom);
        /* TODO: can be optimized as the ring would have the
         *       same vertices, possibly in different order.
         *       maybe comparing number of points could already be
         *       useful.
         */
        if ( GEOSEquals_r(handle, f2er, hole) ) {
#ifdef RTDEBUGF
          RTDEBUGF(ctx, 2, "Hole %d/%d of face %d is face %d", h+1, nholes, i, j);
#endif
          f2->parent = f;
          break;
        }
      }
    }
  }
}

static GEOSGeometry*
collectFacesWithEvenAncestors(GEOSContextHandle_t handle, Face** faces, int nfaces)
{
  GEOSGeometry **geoms = new GEOSGeometry*[nfaces];
  GEOSGeometry *ret;
  unsigned int ngeoms = 0;
  int i;

  for (i=0; i<nfaces; ++i) {
    Face *f = faces[i];
    if ( countParens(handle, f) % 2 ) continue; /* we skip odd parents geoms */
    geoms[ngeoms++] = GEOSGeom_clone_r(handle, f->geom);
  }

  ret = GEOSGeom_createCollection_r(handle, GEOS_MULTIPOLYGON, geoms, ngeoms);
  delete[] geoms;
  return ret;
}

//#define DUMP_GEOM
#ifdef DUMP_GEOM
#include <geos/geom/Geometry.h>
#include <geos/io/WKTWriter.h>

void dumpGeometry(GEOSGeometry* geom)
{
  geos::io::WKTWriter oWriter;
  std::cerr << oWriter.write((geos::geom::Geometry*)geom) << std::endl;
}
#endif

GEOSGeometry*
GEOSBuildArea_r(GEOSContextHandle_t handle, const GEOSGeometry* geom_in)
{
  GEOSGeometry *tmp;
  GEOSGeometry *geos_result, *shp;
  GEOSGeometry const *vgeoms[1];
  uint32_t i, ngeoms;
  int srid = GEOSGetSRID_r(handle, geom_in);
  Face ** geoms;

  vgeoms[0] = geom_in;
#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Polygonizing");
#endif
  geos_result = GEOSPolygonize_r(handle, vgeoms, 1);

#ifdef DUMP_GEOM
  std::cerr << "After polygonize : ";
  dumpGeometry(geos_result);
#endif

#ifdef RTDEBUGF
  RTDEBUGF(ctx, 3, "GEOSpolygonize returned @ %p", geos_result);
#endif
  /* Null return from GEOSpolygonize _r(handle, an exception) */
  if ( ! geos_result ) return 0;

  /*
   * We should now have a collection
   */
#if PARANOIA_LEVEL > 0
  if ( GEOSGeometryTypeId_r(handle, geos_result) != RTCOLLECTIONTYPE )
  {
    GEOSGeom_destroy_r(handle, geos_result);
    rterror(ctx, "Unexpected return from GEOSpolygonize");
    return 0;
  }
#endif

  ngeoms = GEOSGetNumGeometries_r(handle, geos_result);
#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Num geometries from polygonizer: %d", ngeoms);
#endif

#ifdef RTDEBUGF
  RTDEBUGF(ctx, 3, "GEOSpolygonize: ngeoms in polygonize output: %d", ngeoms);
  RTDEBUGF(ctx, 3, "GEOSpolygonize: polygonized:%s",
              rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, geos_result, 0)));
#endif
  /*
   * No geometries in collection, early out
   */
  if ( ngeoms == 0 )
  {
    GEOSSetSRID_r(handle, geos_result, srid);
    return geos_result;
  }

  /*
   * Return first geometry if we only have one in collection,
   * to avoid the unnecessary Geometry clone below.
   */
  if ( ngeoms == 1 )
  {
    tmp = (GEOSGeometry *)GEOSGetGeometryN_r(handle, geos_result, 0);
    if ( ! tmp )
    {
      GEOSGeom_destroy_r(handle, geos_result);
      return 0; /* exception */
    }
    shp = GEOSGeom_clone_r(handle, tmp);
    GEOSGeom_destroy_r(handle, geos_result); /* only safe after the clone above */
    GEOSSetSRID_r(handle, shp, srid);
    return shp;
  }
#ifdef RTDEBUGF
  RTDEBUGF(ctx, 2, "Polygonize returned %d geoms", ngeoms);
#endif
  /*
   * Polygonizer returns a polygon for each face in the built topology.
   *
   * This means that for any face with holes we'll have other faces
   * representing each hole. We can imagine a parent-child relationship
   * between these faces.
   *
   * In order to maximize the number of visible rings in output we
   * only use those faces which have an even number of parents.
   *
   * Example:
   *
   *   +---------------+
   *   |     L0        |  L0 has no parents
   *   |  +---------+  |
   *   |  |   L1    |  |  L1 is an hole of L0
   *   |  |  +---+  |  |
   *   |  |  |L2 |  |  |  L2 is an hole of L1 (which is an hole of L0)
   *   |  |  |   |  |  |
   *   |  |  +---+  |  |
   *   |  +---------+  |
   *   |               |
   *   +---------------+
   *
   * See http://trac.osgeo.org/postgis/ticket/1806
   *
   */

#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Preparing face structures");
#endif

  /* Prepare face structures for later analysis */
  geoms = new Face*[ngeoms];
  for (i=0; i<ngeoms; ++i)
    geoms[i] = newFace(handle, GEOSGetGeometryN_r(handle, geos_result, i));

#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Finding face holes");
#endif

  /* Find faces representing other faces holes */
  findFaceHoles(handle, geoms, ngeoms);

#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Colletting even ancestor faces");
#endif

  /* Build a MultiPolygon composed only by faces with an
   * even number of ancestors */
  tmp = collectFacesWithEvenAncestors(handle, geoms, ngeoms);

#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Cleaning up");
#endif

  /* Cleanup face structures */
  for (i=0; i<ngeoms; ++i) delFace(handle, geoms[i]);
  delete[] geoms;

  /* Faces referenced memory owned by geos_result.
   * It is safe to destroy geos_result after deleting them. */
  GEOSGeom_destroy_r(handle, geos_result);

#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Self-unioning");
#endif

#ifdef DUMP_GEOM
  std::cerr << "Before GEOSUnionCascaded_r : ";
  dumpGeometry(tmp);
#endif
  /* Run a single overlay operation to dissolve shared edges */
  shp = GEOSUnionCascaded_r(handle, tmp);
  if ( ! shp )
  {
    GEOSGeom_destroy_r(handle, tmp);
    return 0; /* exception */
  }

#ifdef DUMP_GEOM
  std::cerr << "After GEOSUnionCascaded_r : ";
  dumpGeometry(shp);
#endif
#ifdef RTGEOM_PROFILE_BUILDAREA
  rtnotice(ctx, "Final cleanup");
#endif

  GEOSGeom_destroy_r(handle, tmp);

  GEOSSetSRID_r(handle, shp, srid);

  return shp;
}
