package org.hupo.psi.mi.psicquic.view.webapp.controller.search;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.myfaces.orchestra.conversation.annotations.ConversationName;
import org.apache.myfaces.orchestra.viewController.annotations.PreRenderView;
import org.apache.myfaces.orchestra.viewController.annotations.ViewController;
import org.apache.myfaces.trinidad.context.RequestContext;
import org.apache.myfaces.trinidad.model.SortableModel;
import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.InteractionClusteringService;
import org.hupo.psi.mi.psicquic.clustering.Service;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.dao.DaoException;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.client.PsicquicRegistryClientException;
import org.hupo.psi.mi.psicquic.registry.client.registry.DefaultPsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.registry.client.registry.PsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.view.webapp.application.PsicquicThreadConfig;
import org.hupo.psi.mi.psicquic.view.webapp.controller.BaseController;
import org.hupo.psi.mi.psicquic.view.webapp.controller.clustering.UserJobs;
import org.hupo.psi.mi.psicquic.view.webapp.controller.config.PsicquicViewConfig;
import org.hupo.psi.mi.psicquic.view.webapp.model.ClusteringResultDataModel;
import org.hupo.psi.mi.psicquic.view.webapp.model.PsicquicResultDataModel;
import org.hupo.psi.mi.psicquic.view.webapp.util.UrlUtils;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import javax.annotation.PostConstruct;
import javax.faces.context.FacesContext;
import javax.faces.event.ActionEvent;
import javax.servlet.http.HttpServletRequest;
import java.io.IOException;
import java.net.Proxy;
import java.net.SocketTimeoutException;
import java.util.*;
import java.util.concurrent.*;

/**
 * Search controller.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */

@Controller("searchBean")
@Scope("conversation.access")
@ConversationName("general")
@ViewController(viewIds = {"/main.xhtml"})
public class SearchController extends BaseController {

	private static final Log log = LogFactory.getLog(SearchController.class);
	private static final int MAX_NETWORK_SIZE = 50000;
	private static final String CLUSTER_JOB_ID = "clusterJobId";
	private static final String QUERY = "query";
	private static final String INCLUDED = "included";
	private static final String EXCLUDED = "excluded";

	private static final int UNEXPECTED_ERROR = -2;
	private static final int TIME_OUT_EXCEPTION = -1;

	private int threadTimeOut = 5;

	@Autowired
	private transient ApplicationContext applicationContext;

	@Autowired
	private UserQuery userQuery;

	@Autowired
	private PsicquicViewConfig config;

	@Autowired
	private ClusteringContext clusteringContext;

	@Autowired
	private InteractionClusteringService clusteringService;

	@Autowired
	private UserJobs userJobs;

	private List<ServiceType> services;
	private List<ServiceType> allServices;
	private Map<String, ServiceType> servicesMap;
	private Map<String, Boolean> serviceSelectionMap;

	private int totalResults = -1;

	private Map<String, String> activeServices;
	private Map<String, String> inactiveServices;
	private Map<String, SortableModel> resultDataModelMap;
	private Map<String, Integer> resultCountMap;

	private Map<String, List<String>> servicesFormatsMap;
	private Map<String, String> servicesLatestMitabVersion;

	private String[] includedServices;
	private String[] excludedServices;

	private String selectedServiceName;

	private boolean clusterSelected;
	private ClusteringJob job;

	private String mitabUrl;

	private Date registryTimestamp;

	private List<Future> runningTasks;

	private static final int XML_INTERACTIONS_LIMIT = 500;



	public SearchController() {
		this.clusterSelected = false;

		this.serviceSelectionMap = Collections.synchronizedMap(new HashMap<String, Boolean>());
		this.resultDataModelMap = Collections.synchronizedMap(new HashMap<String, SortableModel>());
		this.resultCountMap = Collections.synchronizedMap(new HashMap<String, Integer>());
		this.activeServices = Collections.synchronizedMap(new HashMap<String, String>());
		this.servicesFormatsMap = Collections.synchronizedMap(new HashMap<String, List<String>>());
		this.servicesLatestMitabVersion = Collections.synchronizedMap(new HashMap<String, String>());
		this.inactiveServices = Collections.synchronizedMap(new HashMap<String, String>());
	}

	@PostConstruct
	public void refresh() {
		doPsicquicBinarySearch("*");
	}

	private boolean isPartialRequest() {
		FacesContext context = FacesContext.getCurrentInstance();
		return RequestContext.getCurrentInstance().isPartialRequest(context);
	}

	@PreRenderView
	public void preRender() {

		if (isPartialRequest()) {
			return;
		}

		// TODO add new param 'clusterJobId' that can be combined with parameter 'query'

		FacesContext context = FacesContext.getCurrentInstance();

		String statusParam = context.getExternalContext().getRequestParameterMap().get("status");

		if (statusParam != null && "exp".equals(statusParam)) {
			addInfoMessage("Session Expired", "The session have been restored");
		}

		String jobId = context.getExternalContext().getRequestParameterMap().get(CLUSTER_JOB_ID);
		if (jobId != null) {
			final JobDao jobDao = clusteringContext.getDaoFactory().getJobDao();
			final ClusteringJob job;
			try {
				job = jobDao.getJob(jobId);
				if (job == null) {
					log.error("Could not find job by id: " + jobId);
				} else {
					this.userJobs.getCurrentJobs().add(job);
					this.job = job;
					this.clusterSelected = true;
				}
			} catch (DaoException e) {
				final String msg = "Failed to retrieve requested clustered query";
				addErrorMessage("Error", msg);
				log.error(msg + ": " + jobId, e);
			}
		}

		final String queryParam = context.getExternalContext().getRequestParameterMap().get(QUERY);
		if (queryParam != null && queryParam.length() > 0) {
			userQuery.reset();
			userQuery.setSearchQuery(queryParam);
		}

		services = Collections.synchronizedList(new ArrayList<ServiceType>(allServices));

		// included/excluded services
		final String includedServicesParam = context.getExternalContext().getRequestParameterMap().get(INCLUDED);

		if (includedServicesParam != null && includedServicesParam.length() > 0) {
			processIncludedServices(includedServicesParam);
		}

		final String excludedServicesParam = context.getExternalContext().getRequestParameterMap().get(EXCLUDED);

		if (excludedServicesParam != null && excludedServicesParam.length() > 0) {
			processExcludedServices(excludedServicesParam);
		}

		// search
		if (queryParam != null || includedServicesParam != null || excludedServices != null) {
			doNewBinarySearch();
		}
	}

	public void refreshApp(ActionEvent evt) {
		PsicquicRegistryClient registryClient = new DefaultPsicquicRegistryClient(config.getRegistryURL());
		try {
			Date newRegistryTimeStamp = registryClient.registryTimestamp();

			if (registryTimestamp != null && (newRegistryTimeStamp.equals(registryTimestamp) || newRegistryTimeStamp.before(registryTimestamp))) {
				return;
			} else {
				registryTimestamp = newRegistryTimeStamp;
			}
		} catch (PsicquicRegistryClientException e) {
			log.error("Problem checking timestamp from the Registry", e);
			addErrorMessage("Problem checking timestamp from the Registry","Please try again later");
//			throw new RuntimeException("Problem checking timestamp from the Registry", e);
		}

		this.resultDataModelMap.clear();
		this.resultCountMap.clear();
		this.activeServices.clear();
		this.servicesFormatsMap.clear();
		this.servicesLatestMitabVersion.clear();
		this.inactiveServices.clear();

		try {
			refreshServices();
		} catch (PsicquicRegistryClientException e) {
			log.error("Problem checking timestamp from the Registry", e);
			addErrorMessage("Problem checking timestamp from the Registry","Please try again later");
//			throw new RuntimeException("Problem loading services from the Registry", e);
		}
	}

	public String doNewBinarySearch() {

		String searchQuery = userQuery.getSearchQuery();

		if (!clusterSelected) {
			doPsicquicBinarySearch(searchQuery);
		} else {
			doClusteredJobBinarySearch(searchQuery);
		}

		return "interactions";
	}

	public String getActiveServicesName() {
		StringBuilder sb = new StringBuilder(256);
		final Iterator<String> it = activeServices.keySet().iterator();
		while (it.hasNext()) {
			String serviceName = it.next();
			sb.append(serviceName);
			if (it.hasNext()) {
				sb.append(", ");
			}
		}
		return sb.toString();
	}

	public String getInactiveServicesName() {
		StringBuilder sb = new StringBuilder(256);
		final Iterator<String> it = inactiveServices.keySet().iterator();
		while (it.hasNext()) {
			String serviceName = it.next();
			sb.append(serviceName);
			if (it.hasNext()) {
				sb.append(", ");
			}
		}
		return sb.toString();
	}

	private void doClusteredJobBinarySearch(String searchQuery) {
		if (log.isDebugEnabled()) {
			log.debug("\tcluster query:  " + searchQuery);
		}

		final String jobId = job.getJobId();
		final String serviceName = "Clustered query: '" + job.getMiql() + "' from " + printServiceNames(job.getServices());
		selectedServiceName = serviceName;

		final ServiceType service = new ServiceType();
		service.setName(serviceName);

		// Enable download of clustered data via table via a servlet
		service.setRestUrl("/clusterDownload?jobId=" + jobId + "&query=*");

		servicesMap.put(serviceName, service);
		int totalCount = job.getClusteredInteractionCount();
		resultCountMap.put(serviceName, totalCount);
		totalResults = totalCount;

		loadResults(service);
	}

	private String printServiceNames(List<Service> services) {
		StringBuilder sb = new StringBuilder(256);
		final Iterator<Service> it = services.iterator();
		while (it.hasNext()) {
			Service service = it.next();
			sb.append(service.getName());
			if (it.hasNext()) {
				sb.append(", ");
			}
		}
		return sb.toString();
	}

	public void doPsicquicBinarySearch(String searchQuery) {

		if (searchQuery == null) {
			log.warn("The query is null");
			return;
		}

		if (log.isDebugEnabled()) {
			log.debug("\tquery:  " + searchQuery);
		}

		// A minimum checking of the query
		if (searchQuery.length() > 1 && (searchQuery.trim().startsWith("*") || searchQuery.startsWith("?"))) {
			userQuery.setSearchQuery("*:*");
			final String msg = "Your query '" + searchQuery + "' is not correctly formatted";
			final String detail = "Currently we do not support queries prefixed with wildcard characters such as '*' or '?'. " +
					"However, wildcard characters can be used anywhere else in a query (eg. g?vin or gav* for gavin). " +
					"Please do reformat your query.";
			addErrorMessage(msg, detail);
			log.error(msg);
		} else {
			if (isIncludeNegativeQuery()) {
				final String msg = "Your results may contain negative interactions.";
				final String detail = "These results will appear in red rows and they will be excluded from the the " +
						"cluster algorithm and the graph. They can only be exported in MITAB 2.6, MITAB 2.7 and XML 2.5";
				addInfoMessage(msg, detail);
				log.info(msg);
			}
			// searchAndCreateResultModels calls:
			// refreshApp that clean all the maps
			// countInPsicquicService
			searchAndCreateResultModels();

			totalResults = 0;

			for (int count : resultCountMap.values()) {
				if (count > 0) {
					totalResults += count;
				}
			}

			if (totalResults == 0) {
				final String msg = "Your query didn't return any results";
				final String detail = "Use a different query";
				addWarningMessage(msg ,detail);
				log.warn(msg);
			}

			setSelectedServiceName(null);
		}
	}

	public void loadResults(ServiceType service) {
		generateMitabUrl();

		SortableModel results = null;
		try {
			if (!clusterSelected) {
				PsicquicSimpleClient psicquicSimpleClient;
				Proxy proxy = config.getProxy();
				if(proxy != null){
					psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl(), proxy);
				} else {
					psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl());
				}

				results = new PsicquicResultDataModel(
						psicquicSimpleClient,
						getUserQuery().getFilteredSearchQuery(),
						servicesLatestMitabVersion.get(service.getName()));
			} else {
				// TODO for the time being we load ALL clustered data, later we could allow further filtering
					results = new ClusteringResultDataModel(clusteringService, job, "*");
			}
			resultDataModelMap.put(service.getName(), results);
		} catch (IOException e) {
			log.error("Error while building results", e);
		}
	}


	private void generateMitabUrl() {
		final String serviceName = getSelectedServiceName();
		final ServiceType serviceType = getServicesMap().get(serviceName);
		final String query = UrlUtils.encodeUrl(getUserQuery().getFilteredSearchQuery());
		final String restUrl = serviceType.getRestUrl();
		final String queryUrl;

		if (!clusterSelected) {
			boolean endWithSlash = restUrl.endsWith("/");
			queryUrl = restUrl + (endWithSlash ? "" : "/") + "query/" + query;
		} else {
			queryUrl = getApplicationUrl() + restUrl;
		}

		if (log.isDebugEnabled()) log.debug("Reading data from: " + queryUrl);

		mitabUrl = queryUrl;
	}

	private void refreshServices() throws PsicquicRegistryClientException {

		if (allServices == null) {
			PsicquicRegistryClient registryClient = new DefaultPsicquicRegistryClient(config.getRegistryURL());

			if (config.getRegistryTagsAsString() != null) {
				allServices = registryClient.listServices("STATUS", true, config.getRegistryTagsAsString());
			} else {
				allServices = registryClient.listServices();
			}
		}

		services = new ArrayList<ServiceType>(allServices);

		String included = config.getIncludedServices();
		String excluded = config.getExcludedServices();

		if (included != null && included.length() > 0) {
			processIncludedServices(included);
		} else if (excluded != null && excluded.length() > 0) {
			processExcludedServices(excluded);
		}

		populateActiveServicesMap();
		populateInactiveServicesMap();
		populateActiveServicesFormatsAndVersionMap();

		servicesMap = new HashMap<String, ServiceType>();

		boolean populateSelectionMap = serviceSelectionMap.isEmpty();

		for (ServiceType service : services) {
			servicesMap.put(service.getName(), service);

			if (populateSelectionMap) {
				serviceSelectionMap.put(service.getName(), service.isActive());
			}
		}
	}

	private void searchAndCreateResultModels() {
		refreshApp(null);

		final String filteredSearchQuery = userQuery.getFilteredSearchQuery();

		// count the results
		PsicquicThreadConfig threadConfig = (PsicquicThreadConfig) applicationContext.getBean("psicquicThreadConfig");

		ExecutorService executorService = threadConfig.getExecutorService();

		if (runningTasks == null) {
			runningTasks = new ArrayList<Future>();
		} else {
			runningTasks.clear();
		}

		for (final ServiceType service : services) {
			if (service.isActive() && serviceSelectionMap.get(service.getName()) != null && serviceSelectionMap.get(service.getName())) {
				Runnable runnable = new Runnable() {
					public void run() {
						int count = countInPsicquicService(service, filteredSearchQuery);
						resultCountMap.put(service.getName(), count);
					}
				};
				runningTasks.add(executorService.submit(runnable));
			}
		}

		checkAndResumePsicquicTasks();
	}

	private String collectSelectedServices() {
		StringBuilder buffer = new StringBuilder(1064);

		for (Map.Entry<String, Boolean> entry : serviceSelectionMap.entrySet()) {
			if (entry.getValue()) {
				buffer.append(entry.getKey());
			}
		}
		return buffer.toString();
	}

	private void checkAndResumePsicquicTasks() {

		for (Future f : runningTasks) {
			try {
				f.get(threadTimeOut, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				log.error("The psicquic task was interrupted, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (ExecutionException e) {
				log.error("The psicquic task could not be executed, we cancel the task.", e);
				if (!f.isCancelled()) {
					f.cancel(true);
				}
			} catch (TimeoutException e) {
				log.error("Service task stopped because of time out " + threadTimeOut + "seconds.", e);

				if (!f.isCancelled()) {
					f.cancel(true);
				}
			}
		}

		runningTasks.clear();
	}

	private int countInPsicquicService(ServiceType service, String query) {

		int psicquicCount = 0;

		PsicquicSimpleClient psicquicSimpleClient;

		Proxy proxy = config.getProxy();
		if(proxy != null){
			psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl(), proxy);
		} else {
			psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl());
		}
		psicquicSimpleClient.setReadTimeout(5000);


		try {
			psicquicCount = Long.valueOf(psicquicSimpleClient.countByQuery(query)).intValue();
		} catch (IOException e) {
			if(e instanceof SocketTimeoutException){
				psicquicCount = TIME_OUT_EXCEPTION;
			}
			else {
				psicquicCount = UNEXPECTED_ERROR;
			}
			//TODO Check if it is a 2.7 query and explain
			log.warn("Error while building results based on the query: '"
					+ query + "'" + " SERVICE ("+ service.getName() +") RESPONSE: " + e.getMessage());

		}

		return psicquicCount;

	}

	private void populateActiveServicesFormatsAndVersionMap() {
		for (ServiceType service : services) {
			if (service.isActive()) {

				PsicquicSimpleClient psicquicSimpleClient;
				Proxy proxy = config.getProxy();
				if(proxy != null){
					psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl(), proxy);
				} else {
					psicquicSimpleClient = new PsicquicSimpleClient(service.getRestUrl());
				}
				psicquicSimpleClient.setReadTimeout(5000);

				try {
					List<String> formats = psicquicSimpleClient.getFormats();
					servicesFormatsMap.put(service.getName(), formats);
					servicesLatestMitabVersion.put(service.getName(), getLatestMitabVersion(formats));
				} catch (IOException e) {
					//TODO Add a different color for this services and provide a message
					log.error(service.getName() + ": " + e.getMessage());
					service.setActive(false);
					activeServices.remove(service.getName());
					inactiveServices.put(service.getName(), service.getRestUrl());
					serviceSelectionMap.put(service.getName(), false);
				}
			}
		}
	}

	//TODO Improve the way of obtaing the latest version, it is not too flexible.
	private String getLatestMitabVersion(List<String> formats) {

		//By default the version is tab25
		String sVersion = PsicquicSimpleClient.MITAB25;
		int version = 25;

		for (String format : formats) {
			if (format.contentEquals(PsicquicSimpleClient.MITAB26) && (26 > version)) {
				sVersion = PsicquicSimpleClient.MITAB26;
				version = 26;
			} else if (format.contentEquals(PsicquicSimpleClient.MITAB27) && (27 > version)) {
				sVersion = PsicquicSimpleClient.MITAB27;
				version = 27;
			}
		}

		return sVersion;
	}

	private void populateActiveServicesMap() {
		activeServices.clear();
		for (ServiceType service : services) {
			if (service.isActive()) {
				activeServices.put(service.getName(), service.getRestUrl());
//                activeServices.put(service.getName(), service.getSoapUrl());
			}
		}
	}

	private void populateInactiveServicesMap() {
		inactiveServices.clear();
		for (ServiceType service : services) {
			if (!service.isActive()) {
				inactiveServices.put(service.getName(), service.getRestUrl());
//                inactiveServices.put(service.getName(), service.getSoapUrl());
			}
		}
	}

	private void processExcludedServices(String excludedServicesParam) {
		excludedServices = excludedServicesParam.split(",");

		List<ServiceType> includedServicesList = new ArrayList<ServiceType>(excludedServices.length);

		for (ServiceType service : allServices) {
			boolean excluded = false;
			for (String serviceName : excludedServices) {
				if (serviceName.trim().equalsIgnoreCase(service.getName())) {
					excluded = true;
					break;
				}
			}

			if (!excluded) includedServicesList.add(service);
		}

		services = includedServicesList;
	}

	private void processIncludedServices(String includedServicesParam) {
		includedServices = includedServicesParam.split(",");

		List<ServiceType> includedServicesList = new ArrayList<ServiceType>(includedServices.length);

		for (ServiceType service : allServices) {
			for (String serviceName : includedServices) {
				if (serviceName.trim().equalsIgnoreCase(service.getName())) {
					includedServicesList.add(service);
				}
			}
		}

		services = includedServicesList;
	}

	// Getters & Setters
	/////////////////////

	public int getTotalResults() {
		return totalResults;
	}

	public Map<String, Integer> getResultCountMap() {
		return resultCountMap;
	}

	public Map<String, SortableModel> getResultDataModelMap() {
		return resultDataModelMap;
	}

	public Map<String, String> getActiveServices() {
		return activeServices;
	}

	public Map<String, String> getInactiveServices() {
		return inactiveServices;
	}

	public String[] getAllServiceNames() {
		final String[] services = servicesMap.keySet().toArray(new String[servicesMap.size()]);
		Arrays.sort(services, new ServiceNameComparator());
		return services;
	}

	public String[] getActiveServiceNames() {
		final String[] services = activeServices.keySet().toArray(new String[activeServices.size()]);
		Arrays.sort(services, new ServiceNameComparator());
		return services;
	}

	public String[] getInactiveServiceNames() {
		final String[] services = inactiveServices.keySet().toArray(new String[inactiveServices.size()]);
		Arrays.sort(services, new ServiceNameComparator());
		return services;
	}

	public boolean isIncludeNegativeQuery() {

		String query = userQuery.getSearchQuery();
		return (!query.isEmpty() && (query.contains("negative:true")
				|| query.contains("negative:(true OR false)")
				|| query.contains("negative:(false OR true)")));
	}


	public boolean isOnlyNegativeQuery() {

		String query = userQuery.getSearchQuery();
		return (!query.isEmpty() && (query.contains("negative:true")));
	}

	public boolean isBigQuery() {
		return (resultCountMap.get(selectedServiceName) > XML_INTERACTIONS_LIMIT);
	}

	private class ServiceNameComparator implements Comparator<String> {
		public int compare(String o1, String o2) {
			return o1.toLowerCase().compareTo(o2.toLowerCase());
		}
	}

	public List<ServiceType> getServices() {
		return services;
	}

	public Map<String, ServiceType> getServicesMap() {
		return servicesMap;
	}

	public String getSelectedServiceName() {
		return selectedServiceName;
	}

	public void setSelectedServiceName(String selectedServiceName) {
		this.selectedServiceName = selectedServiceName;

		if (selectedServiceName != null) {
			loadResults(servicesMap.get(selectedServiceName));
		}
	}

	public boolean isClusterSelected() {
		return clusterSelected;
	}

	public void setClusterSelected(boolean clusterSelected) {
		this.clusterSelected = clusterSelected;
	}

	public ClusteringJob getJob() {
		return job;
	}

	public void setJob(ClusteringJob job) {
		this.job = job;
	}

	public UserJobs getUserJobs() {
		return userJobs;
	}

	public void setUserJobs(UserJobs userJobs) {
		this.userJobs = userJobs;
	}

	public String unselectClusterJob() {
		clusterSelected = false;
		job = null;
		return doNewBinarySearch();
	}

	public UserQuery getUserQuery() {
		return userQuery;
	}

	public void setUserQuery(UserQuery userQuery) {
		this.userQuery = userQuery;
	}

	private String encodeUrl(String s) {
		return s.replaceAll(" ", "%20");
	}

	public String getMitabUrl() {
		return mitabUrl;
	}

	public String getApplicationUrl() {

		if (isPartialRequest()) {
			return null;
		}

		final FacesContext facesContext = FacesContext.getCurrentInstance();
		HttpServletRequest request = (HttpServletRequest) facesContext.getExternalContext().getRequest();

		final String scheme = request.getScheme();
		final String serverName = request.getServerName();
		final int serverport = request.getServerPort();
		final String contextPath = request.getContextPath();

		return scheme + "://" + serverName + ":" + serverport + contextPath;
	}

	public boolean isGraphEnabled() {
		Boolean enabled = false;
		final String serviceName = getSelectedServiceName();
		if (serviceName != null) {
			final Integer count = getResultCountMap().get(serviceName);
			if (count != null) {
				enabled = (count <= MAX_NETWORK_SIZE);
			}
		}
		return enabled;
	}

	public int getMaxNetworkSize() {
		return MAX_NETWORK_SIZE;
	}

	public int getSelectedServicesCount() {
		int count = 0;

		for (boolean selected : serviceSelectionMap.values()) {
			if (selected) count++;
		}

		return count;
	}

	public Map<String, Boolean> getServiceSelectionMap() {
		return serviceSelectionMap;
	}

	public void setServiceSelectionMap(Map<String, Boolean> serviceSelectionMap) {
		this.serviceSelectionMap = serviceSelectionMap;
	}

	public int getThreadTimeOut() {
		return threadTimeOut;
	}

	public void setThreadTimeOut(int threadTimeOut) {
		this.threadTimeOut = threadTimeOut;
	}

	public Map<String, List<String>> getFormatsMap() {
		return this.servicesFormatsMap;
	}

	public void setServicesFormatsMap(Map<String, List<String>> servicesFormatsMap) {
		this.servicesFormatsMap = servicesFormatsMap;
	}

	public Integer getXmlLimit() {
		return XML_INTERACTIONS_LIMIT;
	}

	public ApplicationContext getApplicationContext() {
		return applicationContext;
	}

	public void setApplicationContext(ApplicationContext applicationContext) {
		this.applicationContext = applicationContext;
	}
}
