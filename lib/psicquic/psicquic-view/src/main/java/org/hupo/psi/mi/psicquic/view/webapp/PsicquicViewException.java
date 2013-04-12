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
package org.hupo.psi.mi.psicquic.view.webapp;

/**
 * Psicquic view exception.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: PsicquicViewException.java 1482 2012-10-18 11:18:29Z noedelta $
 */
public class PsicquicViewException extends RuntimeException
{

    public PsicquicViewException()
    {
        super();
    }

    public PsicquicViewException(String message)
    {
        super(message);
    }

    public PsicquicViewException(String message, Throwable cause)
    {
        super(message, cause);
    }

    public PsicquicViewException(Throwable cause)
    {
        super(cause);
    }
}
