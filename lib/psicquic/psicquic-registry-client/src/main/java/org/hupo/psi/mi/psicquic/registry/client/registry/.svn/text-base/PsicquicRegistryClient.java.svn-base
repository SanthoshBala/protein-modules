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

import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.client.PsicquicRegistryClientException;

import java.util.Date;
import java.util.List;

/**
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public interface PsicquicRegistryClient {

    List<ServiceType> listActiveServices() throws PsicquicRegistryClientException;

    ServiceType getService(String serviceName) throws PsicquicRegistryClientException ;

    List<ServiceType> listInactiveServices() throws PsicquicRegistryClientException ;

    List<ServiceType> listServices() throws PsicquicRegistryClientException ;

    List<ServiceType> listServices(String action, boolean restricted) throws PsicquicRegistryClientException ;

    List<ServiceType> listServices(String action, boolean restricted, String tagExpression) throws PsicquicRegistryClientException ;

    List<ServiceType> listServicesByTags(String tagExpression) throws PsicquicRegistryClientException;

    List<ServiceType> listServices(String action) throws PsicquicRegistryClientException ;

    Date registryTimestamp() throws PsicquicRegistryClientException;
}
