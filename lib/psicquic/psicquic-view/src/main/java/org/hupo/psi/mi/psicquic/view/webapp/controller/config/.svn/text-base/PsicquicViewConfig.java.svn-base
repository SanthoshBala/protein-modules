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
package org.hupo.psi.mi.psicquic.view.webapp.controller.config;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.view.webapp.controller.BaseController;
import org.springframework.beans.factory.InitializingBean;

import javax.faces.event.ActionEvent;
import java.io.*;
import java.net.Proxy;
import java.util.Date;
import java.util.Properties;

/**
 * Application configuration.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class PsicquicViewConfig extends BaseController implements InitializingBean {

	private static final Log log = LogFactory.getLog(PsicquicViewConfig.class);

	public static final String TEMP_DIR = System.getProperty( "java.io.tmpdir" );
	public static final String FILE_SEPARATOR = System.getProperty( "file.separator" );

	public static final int DEFAULT_CLUSTERING_SIZE_LIMIT = 5000;
	public static final String REGISTRY_URL = "http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/";

	private String configFile;
	private String title;
	private String logoUrl;
	private String bannerBackgroundUrl;
	private String registryTagsAsString;
	private String miqlFilterQuery;
	private String includedServices;
	private String excludedServices;
	private String registryURL;
	private int serviceRows;
	private int clusteringSizeLimit;
	private Proxy proxy;
	private String downloadAllLocation;

	public PsicquicViewConfig() {
	}

	public void afterPropertiesSet() throws Exception {
//		if (configFile == null || configFile.length() == 0) {
//			final String tempFile = File.createTempFile("psicquic-view", ".config.properties").toString();
//			log.warn("No config file defined. Created new temporary configuration file: " + tempFile);
//			configFile = tempFile;
//		}
//
//		configFile = configFile.replaceAll("\\\\:", ":");
//
//		File file = new File(configFile);

//		if (!file.exists()) {
//			saveConfigToFile();
//		}

//		loadConfigFromFile(null);

		DetectProxy detectProxy = new DetectProxy();
		proxy = detectProxy.getProxy();
	}

	public void loadConfigFromFile(ActionEvent evt) throws IOException {
		Properties properties = new Properties();
		FileInputStream configInput = new FileInputStream(configFile);
		properties.load(configInput);
		configInput.close();

		title = properties.getProperty("title");
		logoUrl = properties.getProperty("logo.url");
		bannerBackgroundUrl = properties.getProperty("banner.background.url");
		registryTagsAsString = properties.getProperty("registry.tags");
		miqlFilterQuery = properties.getProperty("query.filter");
		includedServices = properties.getProperty("services.included");
		excludedServices = properties.getProperty("services.excluded");

		registryURL = loadProperty(properties, "registry.url", REGISTRY_URL);
		downloadAllLocation = loadProperty(properties, "download.all.storage.location", buildDefaultDir());
		clusteringSizeLimit = loadIntProperty(properties, "clustering.limit.count", DEFAULT_CLUSTERING_SIZE_LIMIT);

		final String strRows = properties.getProperty("services.rows");
		if (strRows != null) serviceRows = Integer.parseInt(strRows);
	}

	private int loadIntProperty(Properties properties, String name, int defaultValue) {
		final String strValue = properties.getProperty(name);
		if (strValue != null) {
			return Integer.parseInt(strValue);
		} else {
			final String sysStrValue = System.getProperty(name);
			if (sysStrValue != null) {
				return Integer.parseInt(sysStrValue);
			} else {
				log.warn("Setting '" + name + "' to default value: " + defaultValue);
				return defaultValue;
			}
		}
	}

	private String loadProperty(Properties properties, String name, String defaultValue) {
		final String strValue = properties.getProperty(name);
		if (strValue != null) {
			return strValue;
		} else {
			final String sysStrValue = System.getProperty(name);
			if (sysStrValue != null) {
				return sysStrValue;
			} else {
				log.warn("Setting '" + name + "' to default value: " + defaultValue);
				return defaultValue;
			}
		}
	}

	public void saveConfigToFile() throws IOException {
		Properties properties = new Properties();

		if (new File(configFile).exists()) {
			properties.load(new FileInputStream(configFile));
		}

		setProperty(properties, "title", title);
		setProperty(properties, "logo.url", logoUrl);
		setProperty(properties, "banner.background.url", bannerBackgroundUrl);
		setProperty(properties, "registry.tags", registryTagsAsString);
		setProperty(properties, "registry.url", registryURL);
		setProperty(properties, "query.filter", miqlFilterQuery);
		setProperty(properties, "services.included", includedServices);
		setProperty(properties, "services.excluded", excludedServices);
		setProperty(properties, "clustering.limit.count", String.valueOf(clusteringSizeLimit));
		setProperty(properties, "services.rows", String.valueOf(serviceRows));
		setProperty(properties, "download.all.storage.location", downloadAllLocation);

		OutputStream outputStream = new FileOutputStream(configFile);
		properties.store(outputStream, "Configuration updated: " + new Date());
		outputStream.close();
	}

	private void setProperty(Properties props, String key, String value) {
		if (value != null) {
			props.put(key, value);
		}
	}

	public void store(ActionEvent evt) {
		try {
			saveConfigToFile();
			super.addInfoMessage("File saved successfully", configFile);
		} catch (IOException e) {
			e.printStackTrace();
			super.addErrorMessage("Problem saving configuration", e.getMessage());
		}
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public String getConfigFile() {
		return configFile;
	}

	public void setConfigFile(String configFile) {
		this.configFile = configFile;
	}

	public String getLogoUrl() {
		return logoUrl;
	}

	public void setLogoUrl(String logoUrl) {
		this.logoUrl = logoUrl;
	}

	public String getRegistryTagsAsString() {
		return registryTagsAsString;
	}

	public void setRegistryTagsAsString(String registryTagsAsString) {
		this.registryTagsAsString = registryTagsAsString;
	}

	public String getMiqlFilterQuery() {
		return miqlFilterQuery;
	}

	public void setMiqlFilterQuery(String miqlFilterQuery) {
		this.miqlFilterQuery = miqlFilterQuery;
	}

	public String getIncludedServices() {
		return includedServices;
	}

	public void setIncludedServices(String includedServices) {
		this.includedServices = includedServices;
	}

	public String getExcludedServices() {
		return excludedServices;
	}

	public void setExcludedServices(String excludedServices) {
		this.excludedServices = excludedServices;
	}

	public String getBannerBackgroundUrl() {
		return bannerBackgroundUrl;
	}

	public void setBannerBackgroundUrl(String bannerBackgroundUrl) {
		this.bannerBackgroundUrl = bannerBackgroundUrl;
	}

	public int getServiceRows() {
		return serviceRows;
	}

	public void setServiceRows(int serviceRows) {
		this.serviceRows = serviceRows;
	}

	public int getClusteringSizeLimit() {
		return clusteringSizeLimit;
	}

	public void setClusteringSizeLimit(int clusteringSizeLimit) {
		this.clusteringSizeLimit = clusteringSizeLimit;
	}

	public String getRegistryURL() {
		return registryURL;
	}

	public void setRegistryURL(String registryURL) {
		this.registryURL = registryURL;
	}

	public Proxy getProxy() {
		return proxy;
	}


	private String buildDefaultDir() {
		String defaultLocation = TEMP_DIR;
		if( ! defaultLocation.endsWith( FILE_SEPARATOR ) ) {
			defaultLocation += FILE_SEPARATOR;
		}
		defaultLocation += "download-cache";
		createDirectoryIfNotExist( defaultLocation );
		return defaultLocation;
	}

	private File createDirectoryIfNotExist( String dir ) {
		File file = new File( dir );
		if(!file.exists()) {
			file.mkdirs();
		}
		return file;
	}

	public String getDownloadAllLocation() {
		if( downloadAllLocation == null ) {
			downloadAllLocation = buildDefaultDir();
		}
		return downloadAllLocation;
	}

	public File getDownloadAllLocationFile() {
		return new File(getDownloadAllLocation());
	}

	public void setDownloadAllLocation( String downloadAllLocation) {
		File file = createDirectoryIfNotExist( downloadAllLocation );
		this.downloadAllLocation = file.getAbsolutePath();
		log.warn( "Setting download all location to: " + this.downloadAllLocation );
	}
}
