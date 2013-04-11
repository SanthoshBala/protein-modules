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
package org.hupo.psi.mi.psicquic.ws.utils;

import psidev.psi.mi.search.Searcher;

/**
 * Class that can be used to index PSI-MI TAB files.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class MitabIndexer {

    public static void main(String[] args) throws Exception {

        if (args.length != 3) {
            throw new IllegalArgumentException("This class needs three arguments: <indexDirectory> <mitabFile> <booleanHasHeader>");
        }

        // directory (index) to be created
        String indexDirectory = args[0];

        // mitab file to index
        String mitabFile = args[1];

        // whether the mitab file has a header or not.
        boolean hasHeader = Boolean.valueOf(args[2]);

        System.out.println("Creating index: "+indexDirectory+" from file: "+mitabFile);

        // creates the index
        Searcher.buildIndex(indexDirectory, mitabFile, true, hasHeader);
    }
}

