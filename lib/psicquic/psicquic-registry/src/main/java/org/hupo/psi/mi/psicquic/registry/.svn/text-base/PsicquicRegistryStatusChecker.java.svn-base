/**
 * Copyright 2011 The European Bioinformatics Institute, and others.
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
package org.hupo.psi.mi.psicquic.registry;

import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.registry.config.PsicquicRegistryConfig;
import org.hupo.psi.mi.psicquic.registry.config.PsicquicRegistryThreadConfig;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.scheduling.annotation.Scheduled;
import org.springframework.stereotype.Controller;

import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.concurrent.*;

/**
 * Checks the state of the services every 5 mins
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
@Controller
public class PsicquicRegistryStatusChecker {

	private static final Log log = LogFactory.getLog(PsicquicRegistryStatusChecker.class);

	@Autowired
	private PsicquicRegistryConfig config;

	@Autowired
	PsicquicRegistryThreadConfig threadConfig;

	private Date lastRefreshed = new Date();

	private int threadTimeOut = 10;

	private List<Future> runningTasks;

	public PsicquicRegistryStatusChecker() {
		this.runningTasks = new ArrayList<Future>(40);
	}

	@Scheduled(fixedDelay = 5 * 60 * 1000) // every 5 mins
	public void refreshServices() {

		ExecutorService executorService = threadConfig.getExecutorService();

		runningTasks.clear();

		for (final ServiceType serviceStatus : config.getRegisteredServices()) {
			Runnable runnable = new Runnable() {
				public void run() {
					checkStatus(serviceStatus);
				}
			};
			runningTasks.add(executorService.submit(runnable));
		}

		checkAndResumeRegistryTasks();

		lastRefreshed = new Date();
	}

	private void checkAndResumeRegistryTasks() {

		for (Future f : runningTasks) {
			try {
				f.get(threadTimeOut, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				log.error("The registry task was interrupted, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (ExecutionException e) {
				log.error("The registry task could not be executed, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (TimeoutException e) {
				log.error("Service task stopped because of time out " + threadTimeOut + " seconds.", e);

				if (!f.isCancelled()) {
					f.cancel(true);
				}
			}
		}

		runningTasks.clear();
	}

	private void checkStatus(ServiceType serviceStatus) {
		try {

			final URL versionUrl = new URL(serviceStatus.getRestUrl() + "version");
			final URL countURL = new URL(serviceStatus.getRestUrl() + "query/*?format=count");

			final HttpURLConnection urlConnection = (HttpURLConnection) versionUrl.openConnection();
			int code = urlConnection.getResponseCode();

			urlConnection.connect();

			if (HttpURLConnection.HTTP_OK == code) {

				serviceStatus.setActive(true);

				final String version;
				final String strCount;

				InputStream contentStream = null;
				InputStream countStream = null;

				try {

					//TODO Add a double check to know if the service is active
					// or not add a catch block for the exceptions

					contentStream = (InputStream) urlConnection.getContent();
					version = IOUtils.toString(contentStream);
					serviceStatus.setVersion(version);

					countStream = countURL.openStream();
					strCount = IOUtils.toString(countStream);
					serviceStatus.setCount(Long.valueOf(strCount));

				} finally {
					if(contentStream!=null){
						contentStream.close();
					}
					if(countStream!=null){
						countStream.close();
					}
				}
			} else {
				serviceStatus.setActive(false);
			}

            urlConnection.disconnect();

		} catch (Throwable e) {
			serviceStatus.setActive(false);
		}


	}

	public Date getLastRefreshed() {
		return lastRefreshed;
	}
}
