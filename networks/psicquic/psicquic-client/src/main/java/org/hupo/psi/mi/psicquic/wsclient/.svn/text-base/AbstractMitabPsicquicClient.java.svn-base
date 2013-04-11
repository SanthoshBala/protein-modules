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

import org.hupo.psi.mi.psicquic.DbRef;
import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.RequestInfo;
import org.hupo.psi.mi.psicquic.ResultInfo;
import org.hupo.psi.mi.psicquic.wsclient.result.MitabSearchResult;
import psidev.psi.mi.tab.model.builder.PsimiTabVersion;

import java.util.List;

/**
 * TODO comment that class header
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public abstract class AbstractMitabPsicquicClient extends AbstractPsicquicClient<MitabSearchResult> {
    
    private static final String RETURN_TYPE_25 = "psi-mi/tab25";
    private static final String RETURN_TYPE_26 = "psi-mi/tab26";
    private static final String RETURN_TYPE_27 = "psi-mi/tab27";

    protected PsimiTabVersion version = PsimiTabVersion.v2_5;

    public AbstractMitabPsicquicClient(String serviceAddress) {
        super(serviceAddress);
    }

    protected AbstractMitabPsicquicClient(String serviceAddress, long timeout) {
        super(serviceAddress, timeout);
    }

    public AbstractMitabPsicquicClient(String serviceAddress, String proxyHost, Integer proxyPort) {
        super(serviceAddress, proxyHost, proxyPort);
    }

    public AbstractMitabPsicquicClient(String serviceAddress, long timeout, String proxyHost, Integer proxyPort) {
        super(serviceAddress, timeout, proxyHost, proxyPort);
    }

    public AbstractMitabPsicquicClient(String serviceAddress, PsimiTabVersion version) {
        super(serviceAddress);
        this.version = version;
    }

    protected AbstractMitabPsicquicClient(String serviceAddress, long timeout, PsimiTabVersion version) {
        super(serviceAddress, timeout);
        this.version = version;
    }

    public AbstractMitabPsicquicClient(String serviceAddress, String proxyHost, Integer proxyPort, PsimiTabVersion version) {
        super(serviceAddress, proxyHost, proxyPort);
        this.version = version;
    }

    public AbstractMitabPsicquicClient(String serviceAddress, long timeout, String proxyHost, Integer proxyPort, PsimiTabVersion version) {
        super(serviceAddress, timeout, proxyHost, proxyPort);
        this.version = version;
    }

    public MitabSearchResult getByQuery(String query, int firstResult, int maxResults) throws PsicquicClientException {
        RequestInfo requestInfo = createRequestInfo(processMitabVersion(), firstResult, maxResults);

        QueryResponse response;
        try {
            response = getService().getByQuery(query, requestInfo);
        } catch (Exception e) {
            throw new PsicquicClientException("There was a problem running the service", e);
        }

        return createSearchResult(response);
    }

    public MitabSearchResult getByInteractor(String identifier, int firstResult, int maxResults) throws PsicquicClientException {
        String RETURN_TYPE = processMitabVersion();
        DbRef dbRef = createDbRef(identifier);
        RequestInfo requestInfo = createRequestInfo(RETURN_TYPE, firstResult, maxResults);

        QueryResponse response;
        try {
            response = getService().getByInteractor(dbRef, requestInfo);
        } catch (Exception e) {
            throw new PsicquicClientException("There was a problem running the service", e);
        }

        return createSearchResult(response);
    }

    public MitabSearchResult getByInteraction(String identifier, int firstResult, int maxResults) throws PsicquicClientException {
        String RETURN_TYPE = processMitabVersion();

        DbRef dbRef = createDbRef(identifier);
        RequestInfo requestInfo = createRequestInfo(RETURN_TYPE, firstResult, maxResults);

        QueryResponse response;
        try {
            response = getService().getByInteraction(dbRef, requestInfo);
        } catch (Exception e) {
            throw new PsicquicClientException("There was a problem running the service", e);
        }

        return createSearchResult(response);
    }

    public MitabSearchResult getByInteractionList(String[] identifiers, int firstResult, int maxResults) throws PsicquicClientException {
        String RETURN_TYPE = processMitabVersion();

        List<DbRef> dbRefs = createDbRefs(identifiers);
        RequestInfo requestInfo = createRequestInfo(RETURN_TYPE, firstResult, maxResults);

        QueryResponse response;
        try {
            response = getService().getByInteractionList(dbRefs, requestInfo);
        } catch (Exception e) {
            throw new PsicquicClientException("There was a problem running the service", e);
        }

        return createSearchResult(response);
    }

    public MitabSearchResult getByInteractorList(String[] identifiers, QueryOperand operand, int firstResult, int maxResults) throws PsicquicClientException {
        String RETURN_TYPE = processMitabVersion();

        List<DbRef> dbRefs = createDbRefs(identifiers);
        RequestInfo requestInfo = createRequestInfo(RETURN_TYPE, firstResult, maxResults);

        QueryResponse response;
        try {
            response = getService().getByInteractorList(dbRefs, requestInfo, operand.toString());
        } catch (Exception e) {
            throw new PsicquicClientException("There was a problem running the service", e);
        }

        return createSearchResult(response);
    }

    private String processMitabVersion() {
        String RETURN_TYPE = RETURN_TYPE_25;

        if (version != null){
            if (version.equals(PsimiTabVersion.v2_6)){
                RETURN_TYPE = RETURN_TYPE_26;
            }
            else if(version.equals(PsimiTabVersion.v2_7)){
                RETURN_TYPE = RETURN_TYPE_27;
            }
        }
        return RETURN_TYPE;
    }

    private MitabSearchResult createSearchResult(QueryResponse response) throws PsicquicClientException {
        String mitab = response.getResultSet().getMitab();

        ResultInfo resultInfo = response.getResultInfo();

        return new MitabSearchResult(mitab, resultInfo.getTotalResults(),
                                                 resultInfo.getFirstResult(),
                                                 resultInfo.getBlockSize());
    }
}
