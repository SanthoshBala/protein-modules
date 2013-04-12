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
package org.hupo.psi.mi.psicquic.wsclient;

/**
 * TODO comment that class header
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: QueryOperand.java 39 2009-04-24 09:40:12Z brunoaranda $
 */
public enum QueryOperand {

    AND("AND"),
    OR("OR");

    private String strOperand;

    private QueryOperand(String strOperand) {
        this.strOperand = strOperand;
    }

    public String toString() {
        return strOperand;
    }
}