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

import psidev.psi.mi.tab.model.CrossReference;

import java.io.IOException;
import java.util.Enumeration;
import java.util.Properties;

/**
 * Xref Link Context Factory.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class XrefLinkContextFactory
{

    private static final String DEFAULT_XREF_LINKS_RES = "/org/hupo/psi/mi/psicquic/binarysearch/default-xref-links.properties";

    private XrefLinkContextFactory() {
    }

    public static XrefLinkContext createDefaultXrefLinkContext() {
        Properties props = new Properties();
        try {
            props.load(XrefLinkContextFactory.class.getResourceAsStream(DEFAULT_XREF_LINKS_RES));
        }
        catch (IOException e) {
            throw new RuntimeException("Problems reading xref-links file from resource: "+DEFAULT_XREF_LINKS_RES);
        }

        XrefLinkContext linkContext = new XrefLinkContext();

        Enumeration<?> propNames = props.propertyNames();

        while (propNames.hasMoreElements()) {
            String propName = (String) propNames.nextElement();
            linkContext.setUrl(propName, props.getProperty(propName));
        }

        return linkContext;
    }

    public static String xrefToString(CrossReference xref, boolean showDb, boolean showId, boolean showText) {
        if (xref == null) return "";

        StringBuilder sb = new StringBuilder( 32 );
        if (showDb) {
            sb.append(xref.getDatabase()).append(":");
        }
        if (showId) {
            sb.append(xref.getIdentifier());
        }
        if (showText && xref.getText() != null) {
            if (showDb || showId) {
              sb.append("(");
            }
            sb.append(xref.getText());

            if (showDb || showId) {
              sb.append(")");
            }
        }
        return sb.toString();
    }
}
