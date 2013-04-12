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
package org.hupo.psi.mi.psicquic.registry.client.registry;

import org.hupo.psi.mi.psicquic.registry.Registry;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.client.PsicquicRegistryClientException;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

/**
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class DefaultPsicquicRegistryClient implements PsicquicRegistryClient {

    private static final String REGISTRY_URL_DEFAULT = "http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/";

    private String registryBaseUrl;

    public DefaultPsicquicRegistryClient() {
        registryBaseUrl = REGISTRY_URL_DEFAULT;
    }

    public DefaultPsicquicRegistryClient(String registryBaseUrl) {
        this.registryBaseUrl = registryBaseUrl;
    }

    public ServiceType getService(String serviceName) throws PsicquicRegistryClientException {
        String url =createDefaultRegistryUrl("STATUS")+"&name="+serviceName;

        try {
            List<ServiceType> serviceTypes = unmarshalUrl(url);

            if (serviceTypes.isEmpty()) {
                return null;
            }

            return serviceTypes.get(0);

        } catch (Throwable t) {
            throw new PsicquicRegistryClientException("Problem reading registry", t);
        }
    }

    public List<ServiceType> listActiveServices() throws PsicquicRegistryClientException {
        return listServices("ACTIVE");
    }

    public List<ServiceType> listServices() throws PsicquicRegistryClientException {
        return listServices("STATUS");
    }

    public List<ServiceType> listInactiveServices() throws PsicquicRegistryClientException {
        return listServices("INACTIVE");
    }

    public List<ServiceType> listServices(String action) throws PsicquicRegistryClientException {
        String url = createDefaultRegistryUrl(action);

        try {
            return unmarshalUrl(url);
        } catch (Throwable t) {
            throw new PsicquicRegistryClientException("Problem reading registry", t);
        }
    }

    public Date registryTimestamp() throws PsicquicRegistryClientException {
        Date timestamp = null;

        try {
            URL url = new URL(registryBaseUrl+"registry/timestamp");

            InputStream is = url.openStream();

            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(is));

            String timestampStr = bufferedReader.readLine();

            timestamp = new SimpleDateFormat("E MMM dd HH:mm:ss z yyyy").parse(timestampStr);
        } catch (Throwable t) {
            throw new PsicquicRegistryClientException("Problem reading timestamp from Registry", t);
        }

        return timestamp;
    }


    public List<ServiceType> listServices(String action, boolean restricted) throws PsicquicRegistryClientException {
        String url = createDefaultRegistryUrl(action)+"&restricted="+(restricted? "y" : "n");

        try {
            return unmarshalUrl(url);
        } catch (Throwable t) {
            throw new PsicquicRegistryClientException("Problem reading registry", t);
        } 
    }

    public List<ServiceType> listServices(String action, boolean restricted, String tagExpression) throws PsicquicRegistryClientException {
        try {
            String url = createDefaultRegistryUrl(action)+"&restricted="+(restricted? "y" : "n")+"&tags="+ URLEncoder.encode(tagExpression, "utf-8");
            return unmarshalUrl(url);
        } catch (Throwable t) {
            throw new PsicquicRegistryClientException("Problem reading registry", t);
        }
    }

    public List<ServiceType> listServicesByTags(String tagExpression) throws PsicquicRegistryClientException {
        return listServices("STATUS", false, tagExpression);
    }

    private String createServiceUrl(String base, String paramUrl) {
        if (!base.endsWith("/")) {
            base = base + "/";
        }
        return base+paramUrl;
    }

    private String createDefaultRegistryUrl(String action) {
        return createServiceUrl(registryBaseUrl, "registry?action="+action.toUpperCase()+"&format=xml");
    }

    private List<ServiceType> unmarshalUrl(String url) throws JAXBException, MalformedURLException {
        final JAXBContext jaxbContext = JAXBContext.newInstance(Registry.class.getPackage().getName());
        final Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
        Registry registry = (Registry) unmarshaller.unmarshal(new URL(url));
        return registry.getServices();
    }

    public static void main(String[] args) throws Exception {
        PsicquicRegistryClient client = new DefaultPsicquicRegistryClient();

        System.out.println("Timestamp: "+client.registryTimestamp());

        List<ServiceType> services = client.listServices();
        
        for (ServiceType service : services) {
            System.out.println(service.getName()+" - "+service.getRestUrl());
        }
    }
}
