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

import javax.faces.component.UIParameter;
import javax.faces.event.ActionEvent;

/**
 * JSF Utilities
 *
 * @author Prem Anand (prem@ebi.ac.uk)
 * @version $Id$
 * @since 2.0.1-SNAPSHOT
 */

public class JsfUtils
{
    private static final String QUOTATION_MARK = "\"";

    /**
     * Uses the converter identified by converterId to convert the value to a String.
     * @value the value to be converted
     * @converterId the id of the converter to be used
     * @componentId the id of the component being rendered
     * @return the String representation of the value.
     */
    public static String valueFromConverter(
            final Object value,
            final String converterId,
            final String componentId)
    {
        final javax.faces.context.FacesContext facesContext=javax.faces.context.FacesContext.getCurrentInstance();
        final javax.faces.convert.Converter converter = facesContext.getApplication().createConverter(converterId);
        return converter.getAsString(facesContext,
                org.apache.commons.lang.StringUtils.isEmpty(componentId)?null:facesContext.getViewRoot().findComponent(componentId),
                value);
    }

    /**
     * Uses the converter identified by converterId to convert the value to a String.
     * @value the value to be converted
     * @converterId the id of the converter to be used
     * @return the String representation of the value.
     */
    public static String valueFromConverter(
            final Object value,
            final String converterId)
    {
        final javax.faces.context.FacesContext facesContext=javax.faces.context.FacesContext.getCurrentInstance();
        final javax.faces.convert.Converter converter = facesContext.getApplication().createConverter(converterId);
        return converter.getAsString(facesContext,null,value);
    }

    /**
     * Returns an javax.faces.event.ActionEvent parameter value, from its name
     * @parameterName the parameter name
     * @event ActionEvent containing the parameter
     * @return the parameter value.
     */
    public static Object getParameterValue(String parameterName, javax.faces.event.ActionEvent event)
    {
        for(Object uiObject : event.getComponent().getChildren())
        {
            if(uiObject instanceof javax.faces.component.UIParameter)
            {
                final javax.faces.component.UIParameter param = (javax.faces.component.UIParameter)uiObject;
                if(param.getName().equals(parameterName))
                {
                    return param.getValue();
                }
            }
        }
        throw new RuntimeException("Parameter "+parameterName+" not found");
    }

    public static String getFirstParamValue(ActionEvent evt) {
        UIParameter param = (UIParameter) evt.getComponent().getChildren().get(0);
        String term = (String) param.getValue();
        return term;
    }

    public static String concat( String a, String b ) {
        return StringUtils.join( new String[]{a, b}, null );
    }

    public static String concat( String a, String b, String c ) {
        return StringUtils.join( new String[]{a, b, c}, null );
    }

    public static String concat( String a, String b, String c, String d ) {
        return StringUtils.join( new String[]{a, b, c, d}, null );
    }

    public static String surroundByQuotesIfMissing(String s) {
        if (s == null) return null;

        if ( ! ( s.startsWith(QUOTATION_MARK) && s.endsWith(QUOTATION_MARK) ) ) {
            s = QUOTATION_MARK + s + QUOTATION_MARK;
        }

        return s;
    }

    /**
     * makes available the current time.
     * @return
     */
    public static long currentTimeMillis() {
        return System.currentTimeMillis();
    }
}
