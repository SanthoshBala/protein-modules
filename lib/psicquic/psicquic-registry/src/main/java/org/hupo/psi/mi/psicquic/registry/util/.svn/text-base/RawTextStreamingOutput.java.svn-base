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
package org.hupo.psi.mi.psicquic.registry.util;

import org.hupo.psi.mi.psicquic.registry.Registry;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

import javax.ws.rs.WebApplicationException;
import javax.ws.rs.core.StreamingOutput;
import java.io.*;

/**
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class RawTextStreamingOutput implements StreamingOutput {

    private Registry registry;

    public RawTextStreamingOutput(Registry registry) {
        this.registry = registry;
    }

    public void write(OutputStream outputStream) throws IOException, WebApplicationException {
        Writer writer = new BufferedWriter(new OutputStreamWriter(outputStream));

        try{
            for (ServiceType service : registry.getServices()) {
                writer.write(service.getName()+"="+service.getSoapUrl()+"\n");
            }
        }
        finally {
            writer.close();
        }
    }
}