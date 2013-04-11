/**
 * Copyright 2007 The European Bioinformatics Institute, and others.
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
 *  limitations under the License.
 */
package org.hupo.psi.mi.psicquic.view.webapp.controller.application;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import javax.faces.context.FacesContext;
import javax.servlet.http.HttpServletRequest;
import java.io.Serializable;

/**
 * Application scope bean, with configuration stuff
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
@Controller
@Scope("request")
public class RequestConfigBean implements Serializable {

    private Log log = LogFactory.getLog(RequestConfigBean.class);

    public static final String DEFAULT_INDEX_LOCATION_INIT_PARAM = "psicquic.DEFAULT_INDEX";

    private String absoluteContextPath;
    private HttpServletRequest request;

    public RequestConfigBean() {
        FacesContext context = FacesContext.getCurrentInstance();
        this.request = (HttpServletRequest) context.getExternalContext().getRequest();

        this.absoluteContextPath = request.getScheme() + "://" +
               request.getServerName() + ":" +
               request.getServerPort() +
               request.getContextPath();
    }

    public HttpServletRequest getRequest() {
        return request;
    }

    public String getAbsoluteContextPath() {
        return absoluteContextPath;
    }
}
