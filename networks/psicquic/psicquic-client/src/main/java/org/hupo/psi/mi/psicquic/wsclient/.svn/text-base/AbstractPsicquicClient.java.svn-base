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

import org.apache.cxf.endpoint.Client;
import org.apache.cxf.frontend.ClientProxy;
import org.apache.cxf.frontend.ClientProxyFactoryBean;
import org.apache.cxf.jaxws.JaxWsProxyFactoryBean;
import org.apache.cxf.transport.common.gzip.GZIPInInterceptor;
import org.apache.cxf.transport.common.gzip.GZIPOutInterceptor;
import org.apache.cxf.transport.http.HTTPConduit;
import org.apache.cxf.transports.http.configuration.HTTPClientPolicy;
import org.hupo.psi.mi.psicquic.DbRef;
import org.hupo.psi.mi.psicquic.PsicquicService;
import org.hupo.psi.mi.psicquic.RequestInfo;

import java.util.ArrayList;
import java.util.List;

/**
 * Abstract superclass for the PSICQUIC clients.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public abstract class AbstractPsicquicClient<T> implements PsicquicClient<T> {

    private PsicquicService service;
    private long connectionTimeOut = 5000L;

    public AbstractPsicquicClient(String serviceAddress) {
        if (serviceAddress == null) return;

        ClientProxyFactoryBean factory = new JaxWsProxyFactoryBean();
        factory.setServiceClass(PsicquicService.class);
        factory.setAddress(serviceAddress);
        factory.getInInterceptors().add(new GZIPInInterceptor());
        factory.getOutInterceptors().add(new GZIPOutInterceptor());

        this.service = (PsicquicService) factory.create();

        final Client client = ClientProxy.getClient(service);

        final HTTPConduit http = (HTTPConduit) client.getConduit();
        final HTTPClientPolicy httpClientPolicy = http.getClient();

        setUpClientProperties(null, null, httpClientPolicy, this.connectionTimeOut);
    }

    public AbstractPsicquicClient(String serviceAddress, long timeout) {
        if (serviceAddress == null) return;

        ClientProxyFactoryBean factory = new JaxWsProxyFactoryBean();
        factory.setServiceClass(PsicquicService.class);
        factory.setAddress(serviceAddress);
        factory.getInInterceptors().add(new GZIPInInterceptor());
        factory.getOutInterceptors().add(new GZIPOutInterceptor());

        this.service = (PsicquicService) factory.create();

        final Client client = ClientProxy.getClient(service);

        final HTTPConduit http = (HTTPConduit) client.getConduit();
        final HTTPClientPolicy httpClientPolicy = http.getClient();

        setUpClientProperties(null, null, httpClientPolicy, timeout);
    }

    public AbstractPsicquicClient(String serviceAddress, String proxyHost, Integer proxyPort) {
        if (serviceAddress == null) return;

        ClientProxyFactoryBean factory = new JaxWsProxyFactoryBean();
        factory.setServiceClass(PsicquicService.class);
        factory.setAddress(serviceAddress);
        factory.getInInterceptors().add(new GZIPInInterceptor());
        factory.getOutInterceptors().add(new GZIPOutInterceptor());

        this.service = (PsicquicService) factory.create();

        final Client client = ClientProxy.getClient(service);

        final HTTPConduit http = (HTTPConduit) client.getConduit();
        final HTTPClientPolicy httpClientPolicy = http.getClient();

        setUpClientProperties(proxyHost, proxyPort, httpClientPolicy, this.connectionTimeOut);
    }

    private void setUpClientProperties(String proxyHost, Integer proxyPort, HTTPClientPolicy httpClientPolicy, long timeout) {
        httpClientPolicy.setReceiveTimeout(timeout);
        httpClientPolicy.setAllowChunking(false);
        httpClientPolicy.setConnectionTimeout(timeout);
        httpClientPolicy.setAcceptEncoding("UTF-8");
        if (proxyHost != null && proxyPort != null){
            httpClientPolicy.setProxyServer(proxyHost);
            httpClientPolicy.setProxyServerPort(proxyPort);
        }
    }

    public AbstractPsicquicClient(String serviceAddress, long timeout, String proxyHost, Integer proxyPort) {
        if (serviceAddress == null) return;

        ClientProxyFactoryBean factory = new JaxWsProxyFactoryBean();
        factory.setServiceClass(PsicquicService.class);
        factory.setAddress(serviceAddress);
        factory.getInInterceptors().add(new GZIPInInterceptor());
        factory.getOutInterceptors().add(new GZIPOutInterceptor());

        this.service = (PsicquicService) factory.create();

        final Client client = ClientProxy.getClient(service);

        final HTTPConduit http = (HTTPConduit) client.getConduit();
        final HTTPClientPolicy httpClientPolicy = http.getClient();

        setUpClientProperties(proxyHost, proxyPort, httpClientPolicy, timeout);
    }

    public PsicquicService getService() {
        return service;
    }

    protected RequestInfo createRequestInfo(String returnType, int firstResult, int maxResults) {
        RequestInfo requestInfo = new RequestInfo();
        requestInfo.setFirstResult(firstResult);
        requestInfo.setBlockSize(maxResults);
        requestInfo.setResultType(returnType);
        return requestInfo;
    }

    protected DbRef createDbRef(String identifier) {
        DbRef dbRef = new DbRef();
        dbRef.setId(identifier);
        return dbRef;
    }

    protected List<DbRef> createDbRefs(String ... identifiers) {
        List<DbRef> dbRefs = new ArrayList<DbRef>(identifiers.length);

        for (String identifier : identifiers) {
            DbRef dbRef = createDbRef(identifier);
            dbRefs.add(dbRef);
        }

        return dbRefs;
    }

    public long getConnectionTimeOut() {
        return connectionTimeOut;
    }

    public void setConnectionTimeOut(long connectionTimeOut) {
        this.connectionTimeOut = connectionTimeOut;
    }
}
