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
package org.hupo.psi.mi.psicquic.view.webapp.util;

import psidev.psi.mi.tab.model.CrossReference;
import org.hupo.psi.mi.psicquic.view.webapp.controller.application.XrefLinkContext;

/**
 * Xref Link Generator.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class XrefLinkGenerator {

    public static String generateUrl(XrefLinkContext linkContext, CrossReference xref) {
        return generateUrl(linkContext, xref, null);
    }

    public static String generateUrl(XrefLinkContext linkContext, CrossReference xref, String linkType) {
        if (linkContext == null) return null;
        if (xref == null) return null;

        return generateUrl(linkContext, xref.getDatabase(), xref.getIdentifier(), linkType);
    }

    public static String generateUrl(XrefLinkContext linkContext, String database, String identifier, String linkType) {
        if (database.equals("taxid")) {
            return replacePlaceholderWithId(linkContext.getOlsUrl(), identifier);
        } else if (database.equals("MI")) {
            return replacePlaceholderWithId(linkContext.getOlsUrl(), "MI:" + identifier);
        } else if (database.equals("pubmed") && !identifier.startsWith("unassigned")) {
            return replacePlaceholderWithId(linkContext.getEuropePubmedCentralUrl(), identifier);
        } else if (database.equals("hierarchView")) {
            return replacePlaceholderWithId(linkContext.getHierarchViewUrl(), identifier);
        } else if ( linkContext.containsKey(database) || linkContext.containsKey(database + "." + linkType) ) {
            // here, the database is the same as the key in the props file (e.g. uniprotkb, intact...)
            String url = null;
            if( linkContext.containsKey( database + "." + linkType ) ){
                url = linkContext.getUrl( database + "." + linkType );
            }
            if( url == null && linkContext.containsKey( database )) {
               url = linkContext.getUrl( database );
            }
            return ( url == null ? null : replacePlaceholderWithId( url, identifier ) );
        }

        return null;
    }

    private static String replacePlaceholderWithId(String url, String identifier) {
        return url.replaceAll("\\{0\\}", identifier);
    }
}
