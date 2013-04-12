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
package org.hupo.psi.mi.psicquic.example;

import org.hupo.psi.mi.psicquic.wsclient.QueryOperand;
import org.hupo.psi.mi.psicquic.wsclient.UniversalPsicquicClient;
import psidev.psi.mi.search.SearchResult;
import psidev.psi.mi.tab.model.BinaryInteraction;

/**
 * Example on how to use the PsicquicClient.
 */
public class PsicquicClientExample {

    public static void main(String[] args) throws Exception {

        // change the enpoint address as needed
        UniversalPsicquicClient client = new UniversalPsicquicClient("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/psicquic");

        SearchResult<BinaryInteraction> searchResult = client.getByQuery("brca2", 0, 200);

        System.out.println("Interactions: "+searchResult.getTotalCount());

    }
}
