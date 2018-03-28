/************************************************************************
 *
 *
 * C-Wrapper for GEOS library
 *
 * Copyright (C) 2010-2012 Sandro Santilli <strk@kbt.io>
 * Copyright (C) 2005-2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 * Author: Sandro Santilli <strk@kbt.io>
 *
 ***********************************************************************/

#ifndef GEOS_UTIL_GEOSCONTEXTHANDLE_H
#define GEOS_UTIL_GEOSCONTEXTHANDLE_H

#include "geos_c.h"

#include <string>

#include <string.h>
#include <stdarg.h>

#include <geos/geom/GeometryFactory.h>
#include <geos/util/Machine.h> // for getMachineByteOrder

using geos::geom::GeometryFactory;

typedef struct GEOSContextHandle_HS
{
    const GeometryFactory *geomFactory;
    char msgBuffer[1024];
    GEOSMessageHandler noticeMessageOld;
    GEOSMessageHandler_r noticeMessageNew;
    void *noticeData;
    GEOSMessageHandler errorMessageOld;
    GEOSMessageHandler_r errorMessageNew;
    void *errorData;
    int WKBOutputDims;
    int WKBByteOrder;
    int initialized;

    GEOSContextHandle_HS()
      :
      geomFactory(0),
      noticeMessageOld(0),
      noticeMessageNew(0),
      noticeData(0),
      errorMessageOld(0),
      errorMessageNew(0),
      errorData(0)
    {
      memset(msgBuffer, 0, sizeof(msgBuffer));
      geomFactory = GeometryFactory::getDefaultInstance();
      WKBOutputDims = 2;
      WKBByteOrder = getMachineByteOrder();
      setNoticeHandler(NULL);
      setErrorHandler(NULL);
      initialized = 1;
    }

    GEOSMessageHandler
    setNoticeHandler(GEOSMessageHandler nf)
    {
        GEOSMessageHandler f = noticeMessageOld;
        noticeMessageOld = nf;
        noticeMessageNew = NULL;
        noticeData = NULL;

        return f;
    }

    GEOSMessageHandler
    setErrorHandler(GEOSMessageHandler nf)
    {
        GEOSMessageHandler f = errorMessageOld;
        errorMessageOld = nf;
        errorMessageNew = NULL;
        errorData = NULL;

        return f;
    }

    GEOSMessageHandler_r
    setNoticeHandler(GEOSMessageHandler_r nf, void *userData) {
        GEOSMessageHandler_r f = noticeMessageNew;
        noticeMessageOld = NULL;
        noticeMessageNew = nf;
        noticeData = userData;

        return f;
    }

    GEOSMessageHandler_r
    setErrorHandler(GEOSMessageHandler_r ef, void *userData)
    {
        GEOSMessageHandler_r f = errorMessageNew;
        errorMessageOld = NULL;
        errorMessageNew = ef;
        errorData = userData;

        return f;
    }

    void
    NOTICE_MESSAGE(const std::string& fmt, ...)
    {
      if (NULL == noticeMessageOld && NULL == noticeMessageNew) {
        return;
      }

      va_list args;
      va_start(args, fmt);
      int result = vsnprintf(msgBuffer, sizeof(msgBuffer) - 1, fmt.c_str(), args);
      va_end(args);

      if (result > 0) {
        if (noticeMessageOld) {
          noticeMessageOld("%s", msgBuffer);
        } else {
          noticeMessageNew(msgBuffer, noticeData);
        }
      }
    }

    void
    ERROR_MESSAGE(const std::string& fmt, ...)
    {
      if (NULL == errorMessageOld && NULL == errorMessageNew) {
        return;
      }

      va_list args;
      va_start(args, fmt);
      int result = vsnprintf(msgBuffer, sizeof(msgBuffer) - 1, fmt.c_str(), args);
      va_end(args);

      if (result > 0) {
        if (errorMessageOld) {
          errorMessageOld("%s", msgBuffer);
        } else {
          errorMessageNew(msgBuffer, errorData);
        }
      }
    }
} GEOSContextHandleInternal_t;

#endif // GEOS_UTIL_GEOSCONTEXTHANDLE_H
