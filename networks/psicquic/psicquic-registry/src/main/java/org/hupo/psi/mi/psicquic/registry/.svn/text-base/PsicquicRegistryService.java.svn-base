/**
 * Copyright 2008 The European Bioinformatics Institute, and others.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hupo.psi.mi.psicquic.registry;

import javax.ws.rs.DefaultValue;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.QueryParam;

/**
 * RESTful web service.
 *
 * v.1.0/search/query/brca2
 * current/search/query/species:human?format=xml25&firstResult=50&maxResults=100
 *
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
@Path("/registry")
public interface PsicquicRegistryService {

    @GET
    @Path("/")
    Object executeAction(@QueryParam("action") String action,
                           @QueryParam("name") String name,
                           @DefaultValue("0") @QueryParam("url") String url,
                           @DefaultValue("html") @QueryParam("format") String format,
                           @DefaultValue("y") @QueryParam("restricted") String showRestricted,
                           @DefaultValue("") @QueryParam("tags") String tags,
                           @DefaultValue("") @QueryParam("excluded") String excludedUrls) throws IllegalActionException;

    @GET
    @Path("/version")
    String getVersion();

    @GET
    @Path("/timestamp")
    String getUpdatedTimestamp();

    @GET
    @Path("/refresh")
    String refresh();

}
