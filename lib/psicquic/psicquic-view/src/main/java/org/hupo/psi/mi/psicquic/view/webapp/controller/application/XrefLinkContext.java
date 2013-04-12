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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * URL Factory.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: XrefLinkContext.java 1777 2013-01-11 11:24:00Z mdumousseau@yahoo.com $
 */
public class XrefLinkContext implements Serializable {

    private Map<String, String> xrefUrls;

    XrefLinkContext() {
        this.xrefUrls = new HashMap<String, String>(16);
    }

    public String getUrl(String xrefKey) {
        xrefKey = prepareKey(xrefKey);

        String url = xrefUrls.get(xrefKey);

        if (url == null) {
            throw new IllegalArgumentException("No URL for key '" + xrefKey + "' found");
        }

        return url;
    }

    public void setUrl(String xrefKey, String url) {
        xrefUrls.put(xrefKey, url);
    }

    public boolean containsKey(String xrefKey) {
        xrefKey = prepareKey(xrefKey);
        return xrefUrls.containsKey(xrefKey);
    }

    private String prepareKey(String xrefKey) {
        xrefKey = xrefKey.replaceAll(" ", "_");
        return xrefKey;
    }

    //
    // Some getters for the most commons cases. It is not necessary to add a getter for each property
    //

    public String getUniprotUrl() {
        return xrefUrls.get("uniprotkb");
    }


    public String getOlsUrl() {
        return xrefUrls.get("ols");
    }

    public String getEuropePubmedCentralUrl() {
        return xrefUrls.get("europePubmedCentral");
    }

    public String getIntactUrl() {
        return xrefUrls.get("intact");
    }

    public String getHierarchViewUrl() {
        return xrefUrls.get("hierarchView");
    }
}
