/**
 * Copyright 2009 The European Bioinformatics Institute, and others.
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
package org.hupo.psi.mi.psicquic.view.webapp.controller.search;

/**
 * A search field.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class SearchField {

    private String name;
    private String displayName;
    private Object autocompleteMethod;

	public SearchField(String name, String displayName) {
        this.name = name;
        this.displayName = displayName;
    }

    public SearchField(String name, String displayName, Object autocompleteMethod) {
        this.name = name;
        this.displayName = displayName;
        this.autocompleteMethod = autocompleteMethod;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getDisplayName() {
        return displayName;
    }

    public void setDisplayName(String displayName) {
        this.displayName = displayName;
    }

    public boolean isAutocompleteCapable() {
        return (autocompleteMethod != null);
    }

    public Object getAutocompleteMethod() {
        return autocompleteMethod;
    }

    public void setAutocompleteMethod(Object autocompleteMethod) {
        this.autocompleteMethod = autocompleteMethod;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SearchField that = (SearchField) o;

        if (displayName != null ? !displayName.equals(that.displayName) : that.displayName != null) return false;
        if (name != null ? !name.equals(that.name) : that.name != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = name != null ? name.hashCode() : 0;
        result = 31 * result + (displayName != null ? displayName.hashCode() : 0);
        return result;
    }
}
