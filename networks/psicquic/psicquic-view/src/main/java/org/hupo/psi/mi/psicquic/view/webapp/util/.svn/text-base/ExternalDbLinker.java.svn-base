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
package org.hupo.psi.mi.psicquic.view.webapp.util;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.myfaces.trinidad.render.ExtendedRenderKitService;
import org.apache.myfaces.trinidad.util.Service;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import javax.faces.context.FacesContext;

/**
 * Utility class for linking to external database resources
 *
 * @author Prem Anand (prem@ebi.ac.uk)
 * @version $Id$
 * @since 2.0.1-SNAPSHOT
 */
@Controller
@Scope( "conversation.access" )
public class ExternalDbLinker {

    public ExternalDbLinker(){

    }

    private static final Log log = LogFactory.getLog( ExternalDbLinker.class );

    //URL Links
    public static final String INTERPROURL = "http://www.ebi.ac.uk/interpro/ISpy?ac=";
    public static final String CHROMOSOMEURL = "http://www.ensembl.org/Homo_sapiens/featureview?type=ProteinAlignFeature;id=";
    public static final String EXPRESSIONURL_PREFIX = "http://www.ebi.ac.uk/microarray-as/atlas/qrs?gprop_0=&gval_0=";
    public static final String EXPRESSIONURL_SUFFIX = "&fexp_0=UP_DOWN&fact_0=&specie_0=&fval_0=(all+conditions)&view=hm";
    public static final String REACTOMEURL = "http://www.reactome.org/cgi-bin/skypainter2";

    //identifier seperators
    public static final String INTERPRO_SEPERATOR = ",";
    public static final String CHROMOSOME_SEPERATOR = ";id=";
    public static final String EXPRESSION_SEPERATOR = "+";
    public static final String REACTOME_SEPERATOR = "\n";


    public void goExternalLink( String baseUrl, String seperator, String[] selected ) {
        goExternalLink( baseUrl, "", seperator, selected );
    }

    public void goExternalLink( String baseUrl, String urlSuffix, String seperator, String[] selected ) {
        FacesContext facesContext = FacesContext.getCurrentInstance();
        ExtendedRenderKitService service = Service.getRenderKitService( facesContext, ExtendedRenderKitService.class );

        if ( selected.length > 0 ) {
            String url = baseUrl + StringUtils.join( selected, seperator ) + urlSuffix;
            service.addScript( facesContext, "window.open('" + url + "');" );
        } else {
            service.addScript( facesContext, "alert('Selection is empty');" );
        }
    }

    //linking to reactome needs a form submit
    public void reactomeLinker( String baseUrl, String separator, String[] selected, String forwardPage ) {
        FacesContext facesContext = FacesContext.getCurrentInstance();
        ExtendedRenderKitService service = Service.getRenderKitService( facesContext, ExtendedRenderKitService.class );

        if ( selected.length > 0 ) {
            service.addScript( facesContext, getReactomeForm( baseUrl,
                                                              StringUtils.join( selected, separator ),
                                                              forwardPage ) );
        } else {
            service.addScript( facesContext, "alert('Selection is empty');" );
        }
    }

    private String getReactomeForm( String baseUrl, String selectedIds, String forwardPage ) {

        StringBuilder sb = new StringBuilder( 1024 );
        sb.append( "reactomeform = document.createElement('form');" );
        sb.append( "reactomeform.method='post';\n" );
        sb.append( "reactomeform.action='" ).append( baseUrl ).append( "';\n" );
        sb.append( "reactomeform.enctype='multipart/form-data';\n" );
        sb.append( "reactomeform.name='skypainter';\n" );
        sb.append( "reactomeform.target='_blank';" );

        createInputHidden(sb, "reactomeform", "QUERY", selectedIds);
        createInputHidden(sb, "reactomeform", "DB", "gk_current");
        createInputHidden(sb, "reactomeform", "SUBMIT", "1");

        sb.append("document.forms[0].parentNode.appendChild(reactomeform);");
        sb.append( "reactomeform.submit();\n" );

        if ( log.isDebugEnabled() ) {
            log.debug("JavaScript to link to  Reactome: \n" +sb.toString() );
        }

        return sb.toString();
    }

    private void createInputHidden( StringBuilder sb, String formName, String name, String value ) {
        sb.append(name).append("Input = document.createElement('input');");
        sb.append(name).append("Input.type = 'hidden';");
        sb.append(name).append("Input.name = '").append(name).append("';");
        sb.append(name).append("Input.value = '").append(value).append("';");
        sb.append(formName).append(".appendChild(");
        sb.append(name).append("Input");
        sb.append(");\n");
    }
}
