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
package org.hupo.psi.mi.psicquic.view.webapp.model;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.myfaces.trinidad.model.SortCriterion;
import org.apache.myfaces.trinidad.model.SortableModel;
import org.hupo.psi.mi.psicquic.view.webapp.util.FlipInteractors;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;
import psidev.psi.mi.search.SearchResult;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.model.BinaryInteraction;
import org.hupo.psi.mi.psicquic.view.webapp.util.FlipInteractors;

import javax.faces.application.FacesMessage;
import javax.faces.context.FacesContext;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.net.SocketTimeoutException;
import java.util.*;

/**
 * DataModel for dar retreived from PSICQUIC services.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class PsicquicResultDataModel extends SortableModel implements Serializable {

	private static final Log log = LogFactory.getLog(PsicquicResultDataModel.class);

	private static String DEFAULT_SORT_COLUMN = "rigid";

	private String query;
	private PsicquicSimpleClient psicquicClient;

	private SearchResult<BinaryInteraction> result;
	private int rowIndex = -1;
	private int firstResult = 0;
	private int maxResults = 50;
	private String version;

	private String sortColumn = DEFAULT_SORT_COLUMN;

	private static final int UNEXPECTED_ERROR = -2;
	private static final int TIME_OUT_EXCEPTION = -1;

	private Boolean ascending = true;

	private Map<String, Boolean> columnSorts;

	private FlipInteractors flipInteractors = new FlipInteractors();

	//TODO Fix the SortableModel
	public PsicquicResultDataModel(PsicquicSimpleClient psicquicClient, String query, String format) throws IOException {
		if (psicquicClient == null) {
			throw new IllegalArgumentException("Trying to create data model with a null Psicquic client");
		}
		if (query == null) {
			throw new IllegalArgumentException("Trying to create data model with a null query");
		}

		this.psicquicClient = psicquicClient;
		this.query = query;
		this.version = format;

		columnSorts = new HashMap<String, Boolean>(16);

		setRowIndex(-1);
		firstResult = 0;
		maxResults = 50;
		fetchResults();

		setWrappedData(result);
	}

	protected void fetchResults() {
		if (query == null) {
			throw new IllegalStateException("Trying to fetch results for a null query");
		}


		if (log.isDebugEnabled()) log.debug("Fetching results: " + query);
		InputStreamReader in = null;
		InputStream inputStream = null;
		int totalCount = -1;
		List<BinaryInteraction> binaryInteractions = new ArrayList<BinaryInteraction>(maxResults);

		try {
			psicquicClient.setReadTimeout(30000);
			totalCount = Long.valueOf(psicquicClient.countByQuery(query)).intValue();

			inputStream = psicquicClient.getByQuery(query, version, firstResult, maxResults);

			in = new InputStreamReader(inputStream);
			PsimiTabReader reader = new PsimiTabReader();
			Iterator<BinaryInteraction> iterator;

			iterator = reader.iterate(in);
			while (iterator.hasNext()) {
				BinaryInteraction interaction = iterator.next();
				flipIfNecessary(interaction, getSearchQuery());
				binaryInteractions.add(interaction);
			}
		} catch (IOException e) {
			binaryInteractions.clear();
			if (e instanceof SocketTimeoutException) {
				totalCount = TIME_OUT_EXCEPTION;
			} else {
				totalCount = UNEXPECTED_ERROR;
			}
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "Temporarily impossible to retrieve results";
			String details = "Error while building results based on the query: '" + query + "'";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}
			log.error(details + " SERVICE RESPONSE: " + e.getMessage());

		} finally {
			try {
				if (in != null) {
					in.close();
				}
				if (inputStream != null) {
					inputStream.close();
				}
			} catch (IOException e) {
				//keep it silence
			}
		}


		result = new SearchResult<BinaryInteraction>(binaryInteractions, totalCount, firstResult, maxResults, null);
	}

	private void flipIfNecessary(BinaryInteraction interaction, String searchQuery) {
		flipInteractors.flipIfNecessary(interaction, searchQuery);
	}

	public int getRowCount() {
		return Long.valueOf(result.getTotalCount()).intValue();
	}

	public Object getRowData() {
		if (result == null) {
			return null;
		}

		if (!isRowWithinResultRange()) {
			firstResult = getRowIndex();

			fetchResults();

		}

		if (!isRowAvailable()) {
//			throw new IllegalArgumentException("row is unavailable");
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "Temporarily impossible to retrieve results";
			String details = "Error while building results based on the query: '" + query + "'";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}
			log.error("Row is unavailable");
			return null;
		}

		final BinaryInteraction binaryInteraction = result.getData().get(rowIndex - firstResult);
		return binaryInteraction;
	}

	public int getRowIndex() {
		return rowIndex;
	}

	public Object getWrappedData() {
		return result;
	}

	public boolean isRowAvailable() {
		if (result == null) {
			return false;
		}

		return rowIndex >= 0 && rowIndex < result.getTotalCount();
	}

	protected boolean isRowWithinResultRange() {
		return (getRowIndex() >= firstResult) && (getRowIndex() < (firstResult + maxResults));
	}

	public void setRowIndex(int rowIndex) {
		if (rowIndex < -1) {
			throw new IllegalArgumentException("illegal rowIndex " + rowIndex);
		}
//		int oldRowIndex = rowIndex;
		this.rowIndex = rowIndex;
//		if (result != null && oldRowIndex != this.rowIndex) {
//			Object data = isRowAvailable() ? getRowData() : null;
//			DataModelEvent event = new DataModelEvent(this, this.rowIndex, data);
//			DataModelListener[] listeners = getDataModelListeners();
//			for (int i = 0; i < listeners.length; i++) {
//				listeners[i].rowSelected(event);
//			}
//		}
	}

	@Override
	public void setSortCriteria(List<SortCriterion> criteria) {
		if ((criteria == null) || (criteria.isEmpty())) {
			this.sortColumn = DEFAULT_SORT_COLUMN;
			columnSorts.clear();
		} else {
			// only use the first criterion
			SortCriterion criterion = criteria.get(0);

			this.sortColumn = criterion.getProperty();

			if (columnSorts.containsKey(sortColumn)) {
				this.ascending = !columnSorts.get(sortColumn);
			} else {
				this.ascending = criterion.isAscending();
			}
			columnSorts.put(sortColumn, ascending);

			if (log.isDebugEnabled())
				log.debug("\tSorting by '" + criterion.getProperty() + "' " + (criterion.isAscending() ? "ASC" : "DESC"));
		}


		fetchResults();

	}


	/**
	 * Checks to see if the underlying collection is sortable by the given property.
	 *
	 * @param property The name of the property to sort the underlying collection by.
	 * @return true, if the property implements java.lang.Comparable
	 */
	@Override
	public boolean isSortable(String property) {
		return true;
	}

	@Override
	public Object getRowKey() {
		return isRowAvailable()
				? getRowIndex()
				: null;
	}

	@Override
	public void setRowKey(Object key) {
		if (key == null) {
			setRowIndex(-1);
		} else {
			setRowIndex((Integer) key);
		}
	}


	public SearchResult<BinaryInteraction> getResult() {
		return result;
	}

	public String getSearchQuery() {
		return query;
	}

	public String getSortColumn() {
		return sortColumn;
	}

	public void setSortColumn(String sortColumn) {
		this.sortColumn = sortColumn;
	}

	public int getFirstResult() {
		return firstResult;
	}

	public boolean isAscending() {
		return ascending;
	}

}
