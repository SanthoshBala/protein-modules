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
package org.hupo.psi.mi.psicquic.registry.config;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.Registry;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.util.RegistryUtil;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

/**
 * Place-holder for the configuration. Initialized by Spring.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class PsicquicRegistryConfig {

    private String version;
    private long serviceCheckTimeout;
    private String registryXmlUrl;
	private List<ServiceType> registeredServices;

	@Autowired
	private ApplicationContext springContext ;


    public PsicquicRegistryConfig() {
        this.serviceCheckTimeout = 1000L;
    }

    @PostConstruct
    public void setup(){
        Registry registry;

        if (registryXmlUrl != null && registryXmlUrl.trim().length() > 0) {
            final URL url;
			InputStream input = null;

			try {
                url = new URL(registryXmlUrl);
				input = url.openStream();
				registry = RegistryUtil.readRegistry(input);

            } catch (Throwable e) {
                throw new RuntimeException("Problem reading registry XML file: "+registryXmlUrl, e);
            } finally {
				if(input!=null){
					try {
						input.close();
					} catch (IOException e) {
						//Keep it quiet
					}
				}
			}
        } else {
            registry = RegistryUtil.readRegistry(PsicquicRegistryConfig.class.getResourceAsStream("/META-INF/psicquic-registry-default.xml"));
        }

        registeredServices = registry.getServices();
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public List<ServiceType> getRegisteredServices() {
        return registeredServices;
    }

    public void setRegisteredServices(List<ServiceType> registeredServices) {
        this.registeredServices = registeredServices;
    }

    public long getServiceCheckTimeout() {
        return serviceCheckTimeout;
    }

    public void setServiceCheckTimeout(long serviceCheckTimeout) {
        this.serviceCheckTimeout = serviceCheckTimeout;
    }

    public String getRegistryXmlUrl() {
		return registryXmlUrl;
	}

	public void setRegistryXmlUrl(String registryXmlUrl) {
		this.registryXmlUrl = registryXmlUrl;
	}

	public void refreshTagOntology(){

		//Cast to the concrete object representation
		SelfDiscoveringOntologyTree miOntologyTree = (SelfDiscoveringOntologyTree)springContext.getBean("miOntologyTree");
		miOntologyTree.clearDiscovered();
	}


	public void refreshRegisteredServices(){

		setup();
	}

}
