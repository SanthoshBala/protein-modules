/**
 * Copyright 2010 The European Bioinformatics Institute, and others.
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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 * This example shows how to access the service and retrieve compressed MITAB results. This uses much less bandwidth, so it may be faster,
 * but increases CPU consumption due to the fact that the response needs to be uncompressed.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class PsicquicSimpleExampleCompression {

    public static void main(String[] args) throws Exception {
        // get a REST URl from the registry http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS
        
        PsicquicSimpleClient client = new PsicquicSimpleClient("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/");

        // count the results
        final long count = client.countByQuery("brca2 OR p53");

        System.out.println("Count: "+count);

        System.out.println("\nMITAB:\n");

        try {
            final InputStream result = client.getByQuery("brca2 OR p53", PsicquicSimpleClient.MITAB25_COMPRESSED);

            // compressed inputstream
            final GZIPInputStream compressedResult = new GZIPInputStream(result);

            BufferedReader in = new BufferedReader(new InputStreamReader(compressedResult));

            String line;

            while ((line = in.readLine()) != null) {
                System.out.println(line);
            }

            in.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
