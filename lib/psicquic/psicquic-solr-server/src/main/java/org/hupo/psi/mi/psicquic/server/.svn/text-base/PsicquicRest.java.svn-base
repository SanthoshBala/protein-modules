package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # RESTful Web service interface
 #
 #=========================================================================== */

import javax.ws.rs.*;

import org.hupo.psi.mi.*;
import org.hupo.psi.mi.psq.*;

public interface PsicquicRest {

    @GET @Path("/interactor/{interactorAc}")
        Object getByInteractor( @PathParam("interactorAc") String intAc,
                                @DefaultValue("") 
                                @QueryParam("db") String db,
                                @DefaultValue("tab25") 
                                @QueryParam("format") String format,
                                @DefaultValue("0") 
                                @QueryParam("firstResult") String firstResult,
                                @DefaultValue("500") 
                                @QueryParam("maxResults") String maxResults) 
        throws PsicquicServiceException,
               NotSupportedMethodException,
               NotSupportedTypeException;
    
    @GET @Path("/interaction/{interactionAc}")
        Object getByInteraction( @PathParam("interactionAc") String intAc,
                                 @DefaultValue("") 
                                 @QueryParam("db") String db,
                                 @DefaultValue("tab25") 
                                 @QueryParam("format") String format,
                                 @DefaultValue("0") 
                                 @QueryParam("firstResult") String firstResult,
                                 @DefaultValue("500") 
                                 @QueryParam("maxResults") String maxResults ) 
        throws PsicquicServiceException,
               NotSupportedMethodException,
               NotSupportedTypeException;

    @GET @Path("/query/{query}")
        Object getByQuery( @PathParam("query") String query,
                           @DefaultValue("tab25") 
                           @QueryParam("format") String format,
                           @DefaultValue("0") 
                           @QueryParam("firstResult") String firstResult,
                           @DefaultValue("500") 
                           @QueryParam("maxResults") String maxResults ) 
        throws PsicquicServiceException,
               NotSupportedMethodException,
               NotSupportedTypeException;

    @POST @Path("/query")
        @Consumes("application/json")
        Object getByPostQuery( String request) 
        throws PsicquicServiceException,
               NotSupportedMethodException,
               NotSupportedTypeException;

    @GET @Path("/formats")
        Object getSupportedFormats() throws PsicquicServiceException,
                                            NotSupportedMethodException;
    
    @GET @Path("/property/{propertyName}")
        Object getProperty( @PathParam("propertyName") 
                            String propertyName )
        throws PsicquicServiceException,
               NotSupportedMethodException;
    
    @GET @Path("/properties")
        Object getProperties() throws PsicquicServiceException,
                                      NotSupportedMethodException;
    
    @GET @Path("/version")
        String getVersion() throws PsicquicServiceException,
                                   NotSupportedMethodException;

}
