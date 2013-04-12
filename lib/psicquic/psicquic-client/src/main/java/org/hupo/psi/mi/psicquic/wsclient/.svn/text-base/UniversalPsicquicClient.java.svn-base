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

import psidev.psi.mi.tab.model.builder.PsimiTabVersion;

/**
 * Universal client for PSICQUIC. It uses the PSI-MI TAB format behind the scenes.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class UniversalPsicquicClient extends AbstractMitabPsicquicClient{

    /**
     * Create a new UniversalPsicquicClient. By default, it does not use a proxy to connect to the webservice and
     * is querying fro mitab 2.5 results
     * @param serviceAddress : address of the webservice
     */
    public UniversalPsicquicClient(String serviceAddress) {
        super(serviceAddress);
    }

    /**
     * Create a new UniversalPsicquicClient. By default, it does not use a proxy to connect to the webservice and
     * is querying fro mitab 2.5 results
     * @param serviceAddress : address of the webservice
     * @param timeout : connection timeout
     */
    public UniversalPsicquicClient(String serviceAddress, long timeout) {
        super(serviceAddress, timeout);
    }

    /**
     * Create a new UniversalPsicquicClient. By default, it is querying fro mitab 2.5 results.
     * The client will use a proxy to connect to the webservice
     * @param serviceAddress : address of the webservice
     * @param proxyHost : proxy host
     * @param proxyPort : proxy port
     */
    public UniversalPsicquicClient(String serviceAddress, String proxyHost, Integer proxyPort) {
        super(serviceAddress, proxyHost, proxyPort);
    }

    /**
     * Create a new UniversalPsicquicClient. By default, it is querying fro mitab 2.5 results.
     * The client will use a proxy to connect to the webservice.
     * @param serviceAddress : address of the webservice
     * @param timeout : connection timeout
     * @param proxyHost : proxy host
     * @param proxyPort : proxy port
     */
    public UniversalPsicquicClient(String serviceAddress, long timeout, String proxyHost, Integer proxyPort) {
        super(serviceAddress, timeout, proxyHost, proxyPort);
    }

    /**
     *
     * @param serviceAddress : address of the webservice
     * @param version : the mitab version of the results (will contain more or less information depending on the service)
     */
    public UniversalPsicquicClient(String serviceAddress, PsimiTabVersion version) {
        super(serviceAddress, version);
    }

    /**
     *
     * @param serviceAddress : address of the webservice
     * @param timeout : connection timeout
     * @param version : the mitab version of the results (will contain more or less information depending on the service)
     */
    protected UniversalPsicquicClient(String serviceAddress, long timeout, PsimiTabVersion version) {
        super(serviceAddress, timeout, version);
    }

    /**
     *
     * @param serviceAddress : address of the webservice
     * @param proxyHost : proxy host
     * @param proxyPort : proxy port
     * @param version : the mitab version of the results (will contain more or less information depending on the service)
     */
    public UniversalPsicquicClient(String serviceAddress, String proxyHost, Integer proxyPort, PsimiTabVersion version) {
        super(serviceAddress, proxyHost, proxyPort, version);
    }

    /**
     *
     * @param serviceAddress : address of the webservice
     * @param timeout : connection timeout
     * @param proxyHost : proxy host
     * @param proxyPort : proxy port
     * @param version : the mitab version of the results (will contain more or less information depending on the service)
     */
    public UniversalPsicquicClient(String serviceAddress, long timeout, String proxyHost, Integer proxyPort, PsimiTabVersion version) {
        super(serviceAddress, timeout, proxyHost, proxyPort, version);
    }
}
