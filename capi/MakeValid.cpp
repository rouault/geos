/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 * 
 * Ported from rtgeom_geos_clean.c from
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
 * Copyright 2009-2010 Sandro Santilli <strk@kbt.io>
 *
 **********************************************************************/

#include "geos_c.h"

#include <geos/util/GEOSContextHandle.h>

/*
 * Return Nth vertex in GEOSGeometry as a POINT.
 * May return NULL if the geometry has NO vertexex.
 */
static GEOSGeometry*
GEOS_getPointN(GEOSContextHandle_t handle,
                      const GEOSGeometry* g_in, uint32_t n)
{
  uint32_t dims;
  const GEOSCoordSequence* seq_in;
  GEOSCoordSeq seq_out;
  double val;
  uint32_t sz;
  int gn;
  GEOSGeometry* ret;

  switch ( GEOSGeomTypeId_r(handle, g_in) )
  {
  case GEOS_MULTIPOINT:
  case GEOS_MULTILINESTRING:
  case GEOS_MULTIPOLYGON:
  case GEOS_GEOMETRYCOLLECTION:
  {
    for (gn=0; gn<GEOSGetNumGeometries_r(handle, g_in); ++gn)
    {
      const GEOSGeometry* g = GEOSGetGeometryN_r(handle, g_in, gn);
      ret = GEOS_getPointN(handle, g,n);
      if ( ret ) return ret;
    }
    break;
  }

  case GEOS_POLYGON:
  {
    ret = GEOS_getPointN(handle, GEOSGetExteriorRing_r(handle, g_in), n);
    if ( ret ) return ret;
    for (gn=0; gn<GEOSGetNumInteriorRings_r(handle, g_in); ++gn)
    {
      const GEOSGeometry* g = GEOSGetInteriorRingN_r(handle, g_in, gn);
      ret = GEOS_getPointN(handle, g, n);
      if ( ret ) return ret;
    }
    break;
  }

  case GEOS_POINT:
  case GEOS_LINESTRING:
  case GEOS_LINEARRING:
    break;

  }

  seq_in = GEOSGeom_getCoordSeq_r(handle, g_in);
  if ( ! seq_in ) return NULL;
  if ( ! GEOSCoordSeq_getSize_r(handle, seq_in, &sz) ) return NULL;
  if ( ! sz ) return NULL;

  if ( ! GEOSCoordSeq_getDimensions_r(handle, seq_in, &dims) ) return NULL;

  seq_out = GEOSCoordSeq_create_r(handle, 1, dims);
  if ( ! seq_out ) return NULL;

  if ( ! GEOSCoordSeq_getX_r(handle, seq_in, n, &val) ) return NULL;
  if ( ! GEOSCoordSeq_setX_r(handle, seq_out, n, val) ) return NULL;
  if ( ! GEOSCoordSeq_getY_r(handle, seq_in, n, &val) ) return NULL;
  if ( ! GEOSCoordSeq_setY_r(handle, seq_out, n, val) ) return NULL;
  if ( dims > 2 )
  {
    if ( ! GEOSCoordSeq_getZ_r(handle, seq_in, n, &val) ) return NULL;
    if ( ! GEOSCoordSeq_setZ_r(handle, seq_out, n, val) ) return NULL;
  }

  return GEOSGeom_createPoint_r(handle, seq_out);
}

/*
 * Fully node given linework
 */
static GEOSGeometry*
GEOS_nodeLines(GEOSContextHandle_t handle, const GEOSGeometry* lines)
{
  GEOSGeometry* noded;
  GEOSGeometry* point;

  /*
   * Union with first geometry point, obtaining full noding
   * and dissolving of duplicated repeated points
   *
   * TODO: substitute this with UnaryUnion?
   */

  point = GEOS_getPointN(handle, lines, 0);
  if ( ! point ) return NULL;

  noded = GEOSUnion_r(handle, lines, point);
  if ( NULL == noded )
  {
    GEOSGeom_destroy_r(handle, point);
    return NULL;
  }

  GEOSGeom_destroy_r(handle, point);

  return noded;
}

/*
 * We expect rtgeom_geos_ensure_init(ctx);
 * Will return NULL on error (expect error handler being called by then)
 *
 */
static GEOSGeometry*
GEOS_makeValidPolygon(GEOSContextHandle_t handle, const GEOSGeometry* gin)
{
  GEOSGeom gout;
  GEOSGeom geos_bound;
  GEOSGeom geos_cut_edges, geos_area, collapse_points;
  GEOSGeometry *vgeoms[3]; /* One for area, one for cut-edges */
  unsigned int nvgeoms=0;

  assert (GEOSGeomTypeId_r(handle, gin) == GEOS_POLYGON ||
          GEOSGeomTypeId_r(handle, gin) == GEOS_MULTIPOLYGON);

  geos_bound = GEOSBoundary_r(handle, gin);
  if ( NULL == geos_bound )
  {
    return NULL;
  }

#ifdef RDDEBUF
  RTDEBUGF(ctx, 3,
                 "Boundaries: %s",
                 rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, geos_bound, 0)));
#endif

  /* Use noded boundaries as initial "cut" edges */

#ifdef RTGEOM_PROFILE_MAKEVALID
  rtnotice(ctx, "ST_MakeValid: noding lines");
#endif


  geos_cut_edges = GEOS_nodeLines(handle, geos_bound);
  if ( NULL == geos_cut_edges )
  {
    GEOSGeom_destroy_r(handle, geos_bound);
    //handle->NOTICE_MESSAGE("RTGEOM_GEOS_nodeLines(ctx): %s", rtgeom_get_last_geos_error(ctx));
    return NULL;
  }

  /* NOTE: the noding process may drop lines collapsing to points.
   *       We want to retrive any of those */
  {
    GEOSGeometry* pi;
    GEOSGeometry* po;

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: extracting unique points from bounds");
#endif

    pi = GEOSGeom_extractUniquePoints_r(handle, geos_bound);
    if ( NULL == pi )
    {
      GEOSGeom_destroy_r(handle, geos_bound);
      //rtnotice(ctx, "GEOSGeom_extractUniquePoints_r(handle): %s",
      //         rtgeom_get_last_geos_error(ctx));
      return NULL;
    }

#ifdef RDDEBUF
    RTDEBUGF(ctx, 3,
                   "Boundaries input points %s",
                   rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, pi, 0)));
#endif

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: extracting unique points from cut_edges");
#endif

    po = GEOSGeom_extractUniquePoints_r(handle, geos_cut_edges);
    if ( NULL == po )
    {
      GEOSGeom_destroy_r(handle, geos_bound);
      GEOSGeom_destroy_r(handle, pi);
      //rtnotice(ctx, "GEOSGeom_extractUniquePoints_r(handle): %s",
      //         rtgeom_get_last_geos_error(ctx));
      return NULL;
    }

#ifdef RDDEBUF
    RTDEBUGF(ctx, 3,
                   "Boundaries output points %s",
                   rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, po, 0)));
#endif

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: find collapse points");
#endif

    collapse_points = GEOSDifference_r(handle, pi, po);
    if ( NULL == collapse_points )
    {
      GEOSGeom_destroy_r(handle, geos_bound);
      GEOSGeom_destroy_r(handle, pi);
      GEOSGeom_destroy_r(handle, po);
      //rtnotice(ctx, "GEOSDifference_r(handle): %s", rtgeom_get_last_geos_error(ctx));
      return NULL;
    }

#ifdef RDDEBUF
    RTDEBUGF(ctx, 3,
                   "Collapse points: %s",
                   rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, collapse_points, 0)));
#endif

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: cleanup(1)");
#endif

    GEOSGeom_destroy_r(handle, pi);
    GEOSGeom_destroy_r(handle, po);
  }
  GEOSGeom_destroy_r(handle, geos_bound);

#ifdef RDDEBUF
  RTDEBUGF(ctx, 3,
                 "Noded Boundaries: %s",
                 rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, geos_cut_edges, 0)));
#endif

  /* And use an empty geometry as initial "area" */
  geos_area = GEOSGeom_createEmptyPolygon_r(handle);
  if ( ! geos_area )
  {
    //rtnotice(ctx, "GEOSGeom_createEmptyPolygon_r(handle): %s", rtgeom_get_last_geos_error(ctx));
    GEOSGeom_destroy_r(handle, geos_cut_edges);
    return NULL;
  }

  /*
   * See if an area can be build with the remaining edges
   * and if it can, symdifference with the original area.
   * Iterate this until no more polygons can be created
   * with left-over edges.
   */
  while (GEOSGetNumGeometries_r(handle, geos_cut_edges))
  {
    GEOSGeometry* new_area=0;
    GEOSGeometry* new_area_bound=0;
    GEOSGeometry* symdif=0;
    GEOSGeometry* new_cut_edges=0;

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: building area from %d edges", GEOSGetNumGeometries_r(handle, geos_cut_edges));
#endif

    /*
     * ASSUMPTION: cut_edges should already be fully noded
     */

    new_area = GEOSBuildArea_r(handle, geos_cut_edges);
    if ( ! new_area )   /* must be an exception */
    {
      GEOSGeom_destroy_r(handle, geos_cut_edges);
      GEOSGeom_destroy_r(handle, geos_area);
      //rtnotice(ctx, "RTGEOM_GEOS_buildArea(ctx) threw an error: %s",
      //         rtgeom_get_last_geos_error(ctx));
      return NULL;
    }

    if ( GEOSisEmpty_r(handle, new_area) )
    {
      /* no more rings can be build with thes edges */
      GEOSGeom_destroy_r(handle, new_area);
      break;
    }

    /*
     * We succeeded in building a ring !
     */

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: ring built with %d cut edges, saving boundaries", GEOSGetNumGeometries_r(handle, geos_cut_edges));
#endif

    /*
     * Save the new ring boundaries first (to compute
     * further cut edges later)
     */
    new_area_bound = GEOSBoundary_r(handle, new_area);
    if ( ! new_area_bound )
    {
      /* We did check for empty area already so
       * this must be some other error */
      //rtnotice(ctx, "GEOSBoundary_r(handle, '%s') threw an error: %s",
      //         rtgeom_to_ewkt(ctx, GEOS2RTGEOM(ctx, new_area, 0)),
      //         rtgeom_get_last_geos_error(ctx));
      GEOSGeom_destroy_r(handle, new_area);
      GEOSGeom_destroy_r(handle, geos_area);
      return NULL;
    }

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: running SymDifference with new area");
#endif

    /*
     * Now symdif new and old area
     */
    symdif = GEOSSymDifference_r(handle, geos_area, new_area);
    if ( ! symdif )   /* must be an exception */
    {
      GEOSGeom_destroy_r(handle, geos_cut_edges);
      GEOSGeom_destroy_r(handle, new_area);
      GEOSGeom_destroy_r(handle, new_area_bound);
      GEOSGeom_destroy_r(handle, geos_area);
      //rtnotice(ctx, "GEOSSymDifference_r(handle) threw an error: %s",
      //         rtgeom_get_last_geos_error(ctx));
      return NULL;
    }

    GEOSGeom_destroy_r(handle, geos_area);
    GEOSGeom_destroy_r(handle, new_area);
    geos_area = symdif;
    symdif = 0;

    /*
     * Now let's re-set geos_cut_edges with what's left
     * from the original boundary.
     * ASSUMPTION: only the previous cut-edges can be
     *             left, so we don't need to reconsider
     *             the whole original boundaries
     *
     * NOTE: this is an expensive operation.
     *
     */

#ifdef RTGEOM_PROFILE_MAKEVALID
    rtnotice(ctx, "ST_MakeValid: computing new cut_edges (GEOSDifference)");
#endif

    new_cut_edges = GEOSDifference_r(handle, geos_cut_edges, new_area_bound);
    GEOSGeom_destroy_r(handle, new_area_bound);
    if ( ! new_cut_edges )   /* an exception ? */
    {
      /* cleanup and throw */
      GEOSGeom_destroy_r(handle, geos_cut_edges);
      GEOSGeom_destroy_r(handle, geos_area);
      /* TODO: Shouldn't this be an rterror ? */
      //rtnotice(ctx, "GEOSDifference_r(handle) threw an error: %s",
      //         rtgeom_get_last_geos_error(ctx));
      return NULL;
    }
    GEOSGeom_destroy_r(handle, geos_cut_edges);
    geos_cut_edges = new_cut_edges;
  }

#ifdef RTGEOM_PROFILE_MAKEVALID
  rtnotice(ctx, "ST_MakeValid: final checks");
#endif

  if ( ! GEOSisEmpty_r(handle, geos_area) )
  {
    vgeoms[nvgeoms++] = geos_area;
  }
  else
  {
    GEOSGeom_destroy_r(handle, geos_area);
  }

  if ( ! GEOSisEmpty_r(handle, geos_cut_edges) )
  {
    vgeoms[nvgeoms++] = geos_cut_edges;
  }
  else
  {
    GEOSGeom_destroy_r(handle, geos_cut_edges);
  }

  if ( ! GEOSisEmpty_r(handle, collapse_points) )
  {
    vgeoms[nvgeoms++] = collapse_points;
  }
  else
  {
    GEOSGeom_destroy_r(handle, collapse_points);
  }

  if ( 1 == nvgeoms )
  {
    /* Return cut edges */
    gout = vgeoms[0];
  }
  else
  {
    /* Collect areas and lines (if any line) */
    gout = GEOSGeom_createCollection_r(handle, GEOS_GEOMETRYCOLLECTION, vgeoms, nvgeoms);
    if ( ! gout )   /* an exception again */
    {
      /* cleanup and throw */
      /* TODO: Shouldn't this be an rterror ? */
      //rtnotice(ctx, "GEOSGeom_createCollection_r(handle) threw an error: %s",
      //         rtgeom_get_last_geos_error(ctx));
      /* TODO: cleanup! */
      return NULL;
    }
  }

  return gout;

}


static GEOSGeometry*
GEOS_makeValidLine(GEOSContextHandle_t handle, const GEOSGeometry* gin)
{
  GEOSGeometry* noded;
  noded = GEOS_nodeLines(handle, gin);
  return noded;
}


static GEOSGeometry*
GEOS_makeValidMultiLine(GEOSContextHandle_t handle, const GEOSGeometry* gin)
{
  GEOSGeometry** lines;
  GEOSGeometry** points;
  GEOSGeometry* mline_out=0;
  GEOSGeometry* mpoint_out=0;
  GEOSGeometry* gout=0;
  uint32_t nlines=0, nlines_alloc;
  uint32_t npoints=0;
  uint32_t ngeoms=0, nsubgeoms;
  uint32_t i, j;

  ngeoms = GEOSGetNumGeometries_r(handle, gin);

  nlines_alloc = ngeoms;
  lines = (GEOSGeometry**)malloc(sizeof(GEOSGeometry*)*nlines_alloc);
  points = (GEOSGeometry**)malloc(sizeof(GEOSGeometry*)*ngeoms);
  if( !lines || !points )
  {
      free(lines);
      free(points);
      return NULL;
  }

  for (i=0; i<ngeoms; ++i)
  {
    const GEOSGeometry* g = GEOSGetGeometryN_r(handle, gin, i);
    GEOSGeometry* vg;
    vg = GEOS_makeValidLine(handle, g);
    if ( GEOSisEmpty_r(handle, vg) )
    {
      /* we don't care about this one */
      GEOSGeom_destroy_r(handle, vg);
    }
    if ( GEOSGeomTypeId_r(handle, vg) == GEOS_POINT )
    {
      points[npoints++] = vg;
    }
    else if ( GEOSGeomTypeId_r(handle, vg) == GEOS_LINESTRING )
    {
      lines[nlines++] = vg;
    }
    else if ( GEOSGeomTypeId_r(handle, vg) == GEOS_MULTILINESTRING )
    {
      nsubgeoms=GEOSGetNumGeometries_r(handle, vg);
      nlines_alloc += nsubgeoms;
      GEOSGeometry** newlines = (GEOSGeometry**)realloc(lines, sizeof(GEOSGeometry*)*nlines_alloc);
      if( !newlines )
      {
          for (i=0; i<nlines; ++i)
          {
              GEOSGeom_destroy_r(handle, lines[i]);
          }
          for (i=0; i<npoints; ++i)
          {
              GEOSGeom_destroy_r(handle, points[i]);
          }
          free(lines);
          free(points);
          return NULL;
      }
      lines = newlines;
      for (j=0; j<nsubgeoms; ++j)
      {
        const GEOSGeometry* gc = GEOSGetGeometryN_r(handle, vg, j);
        /* NOTE: ownership of the cloned geoms will be
         *       taken by final collection */
        lines[nlines++] = GEOSGeom_clone_r(handle, gc);
      }
      GEOSGeom_destroy_r(handle, vg); // was missing in rttopp
    }
    else
    {
      /* NOTE: return from GEOSGeomType will leak
       * but we really don't expect this to happen */
      handle->ERROR_MESSAGE("unexpected geom type returned "
              "by RTGEOM_GEOS_makeValid: %s",
              GEOSGeomType_r(handle, vg));
    }
  }

  if ( npoints )
  {
    if ( npoints > 1 )
    {
      mpoint_out = GEOSGeom_createCollection_r(handle, GEOS_MULTIPOINT,
                                             points, npoints);
    }
    else
    {
      mpoint_out = points[0];
    }
  }

  if ( nlines )
  {
    if ( nlines > 1 )
    {
      mline_out = GEOSGeom_createCollection_r(handle,
                      GEOS_MULTILINESTRING, lines, nlines);
    }
    else
    {
      mline_out = lines[0];
    }
  }

  free(lines);

  if ( mline_out && mpoint_out )
  {
    points[0] = mline_out;
    points[1] = mpoint_out;
    gout = GEOSGeom_createCollection_r(handle, GEOS_GEOMETRYCOLLECTION,
                                     points, 2);
  }
  else if ( mline_out )
  {
    gout = mline_out;
  }
  else if ( mpoint_out )
  {
    gout = mpoint_out;
  }

  free(points);

  return gout;
}

/*
 * We expect rtgeom_geos_ensure_init(ctx);
 * Will return NULL on error (expect error handler being called by then)
 */
static GEOSGeometry*
GEOS_makeValidCollection(GEOSContextHandle_t handle, const GEOSGeometry* gin)
{
  int nvgeoms;
  GEOSGeometry **vgeoms;
  GEOSGeom gout;
  int i;

  nvgeoms = GEOSGetNumGeometries_r(handle, gin);
  if ( nvgeoms == -1 ) {
    //rterror(ctx, "GEOSGetNumGeometries: %s", rtgeom_get_last_geos_error(ctx));
    return 0;
  }

  vgeoms = (GEOSGeometry**)malloc(sizeof(GEOSGeometry*) * nvgeoms );
  if ( ! vgeoms ) {
    handle->ERROR_MESSAGE("RTGEOM_GEOS_makeValidCollection: out of memory");
    return 0;
  }

  for ( i=0; i<nvgeoms; ++i ) {
    vgeoms[i] = GEOSMakeValid_r(handle, GEOSGetGeometryN_r(handle, gin, i) );
    if ( ! vgeoms[i] ) {
      while (i--) GEOSGeom_destroy_r(handle, vgeoms[i]);
      free(vgeoms);
      /* we expect rterror being called already by makeValid */
      return NULL;
    }
  }

  /* Collect areas and lines (if any line) */
  gout = GEOSGeom_createCollection_r(handle, GEOS_GEOMETRYCOLLECTION, vgeoms, nvgeoms);
  if ( ! gout )   /* an exception again */
  {
    /* cleanup and throw */
    for ( i=0; i<nvgeoms; ++i ) GEOSGeom_destroy_r(handle, vgeoms[i]);
    free(vgeoms);
    //rterror(ctx, "GEOSGeom_createCollection_r(handle) threw an error: %s",
    //         rtgeom_get_last_geos_error(ctx));
    return NULL;
  }
  free(vgeoms);

  return gout;

}

GEOSGeometry* GEOSMakeValid_r(GEOSContextHandle_t handle,
                              const GEOSGeometry* gin)
{
    GEOSGeometry* gout;
    char ret_char;

    /*
    * Step 2: return what we got so far if already valid
    */

    ret_char = GEOSisValid_r(handle, gin);
    if ( ret_char == 2 )
    {
        /* I don't think should ever happen */
        return NULL;
    }
    else if ( ret_char )
    {
        /* It's valid at this step, return what we have */
        return GEOSGeom_clone_r(handle, gin);
    }


  /*
   * Step 3 : make what we got valid
   */

    switch (GEOSGeomTypeId_r(handle, gin))
    {
        case GEOS_MULTIPOINT:
        case GEOS_POINT:
            /* points are artays valid, but we might have invalid ordinate values */
            handle->NOTICE_MESSAGE("PUNTUAL geometry resulted invalid to GEOS -- dunno how to clean that up");
            return NULL;
            break;

        case GEOS_LINESTRING:
            gout = GEOS_makeValidLine(handle, gin);
            if ( ! gout )  /* an exception or something */
            {
                return NULL;
            }
            break; /* we've done */

        case GEOS_MULTILINESTRING:
            gout = GEOS_makeValidMultiLine(handle, gin);
            if ( ! gout )  /* an exception or something */
            {
                return NULL;
            }
            break; /* we've done */

        case GEOS_POLYGON:
        case GEOS_MULTIPOLYGON:
        {
            gout = GEOS_makeValidPolygon(handle, gin);
            if ( ! gout )  /* an exception or something */
            {
                return NULL;
            }
            break; /* we've done */
        }

        case GEOS_GEOMETRYCOLLECTION:
        {
            gout = GEOS_makeValidCollection(handle, gin);
            if ( ! gout )  /* an exception or something */
            {
                return NULL;
            }
            break; /* we've done */
        }

        default:
        {
            char* typname = GEOSGeomType_r(handle, gin);
            handle->NOTICE_MESSAGE("ST_MakeValid: doesn't support geometry type: ", typname);
            GEOSFree_r(handle, typname);
            return NULL;
            break;
        }
    }

#if PARANOIA_LEVEL > 1
  /*
   * Now check if every point of input is also found
   * in output, or abort by returning NULL
   *
   * Input geometry was in
   */
  {
      int loss;
      GEOSGeometry *pi, *po, *pd;

      /* TODO: handle some errors here...
       * Lack of exceptions is annoying indeed,
       * I'm getting old --strk;
       */
      pi = GEOSGeom_extractUniquePoints_r(handle, gin);
      po = GEOSGeom_extractUniquePoints_r(handle, gout);
      pd = GEOSDifference_r(handle, pi, po); /* input points - output points */
      GEOSGeom_destroy_r(handle, pi);
      GEOSGeom_destroy_r(handle, po);
      loss = !GEOSisEmpty_r(handle, pd);
      GEOSGeom_destroy_r(handle, pd);
      if ( loss )
      {
        handle->NOTICE_MESSAGE(handle, "Vertices lost in GEOS_makeValid");
        /* return NULL */
      }
  }
#endif /* PARANOIA_LEVEL > 1 */


  return gout;
}
