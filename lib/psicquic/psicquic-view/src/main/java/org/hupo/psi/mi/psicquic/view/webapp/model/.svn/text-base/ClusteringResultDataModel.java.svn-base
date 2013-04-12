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
import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.clustering.ClusteringServiceException;
import org.hupo.psi.mi.psicquic.clustering.DefaultInteractionClusteringService;
import org.hupo.psi.mi.psicquic.clustering.InteractionClusteringService;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.view.webapp.util.FlipInteractors;
import psidev.psi.mi.search.SearchResult;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.model.BinaryInteraction;
import org.hupo.psi.mi.psicquic.view.webapp.util.FlipInteractors;

import javax.faces.application.FacesMessage;
import javax.faces.context.FacesContext;
import javax.faces.model.DataModelEvent;
import javax.faces.model.DataModelListener;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/**
 * DataModel for the clustered data.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 */
public class ClusteringResultDataModel extends SortableModel implements Serializable {

    private static final Log log = LogFactory.getLog( ClusteringResultDataModel.class);

    private static String DEFAULT_SORT_COLUMN = "rigid";

    private ClusteringJob job;

    /**
     * User query on the clustered data.
     */
    private String query;

    private SearchResult<BinaryInteraction> result;
    private int rowIndex = -1;
    private int firstResult = 0;
    private int maxResults = 50;

    private String sortColumn = DEFAULT_SORT_COLUMN;
    private Boolean ascending = true;

    private Map<String,Boolean> columnSorts;

    private InteractionClusteringService clusteringService;

	private FlipInteractors flipInteractors = new FlipInteractors();


	public ClusteringResultDataModel(InteractionClusteringService clusteringService, ClusteringJob job, String query) {
        if (job == null) {
            throw new IllegalArgumentException("Trying to create data model with a null job");
        }
        if (query == null) {
            throw new IllegalArgumentException("Trying to create data model with a null query");
        }
        if ( clusteringService == null ) {
            throw new IllegalArgumentException( "You must give a non null clusteringService" );
        }

        this.clusteringService = clusteringService;
        this.job = job;
        this.query = query;

        columnSorts = new HashMap<String, Boolean>(16);

        setRowIndex(0);
        fetchResults();

        setWrappedData(result);
    }

    protected void fetchResults() {
        if (query == null) {
            throw new IllegalArgumentException("Trying to fetch results for a null query");
        }

        log.debug("Fetching results: " + query);

        // search Lucene index and wrap the data
        final InteractionClusteringService ics = new DefaultInteractionClusteringService();

		//The parameter query in this case retrives all the values in the clustered file
        try {
            final QueryResponse response = clusteringService.query( job.getJobId(),
                                                                    query,
                                                                    firstResult,
                                                                    maxResults,
                                                                    InteractionClusteringService.RETURN_TYPE_MITAB25 );

			String originalQuery = job.getMiql();

            // convert to a list of BinaryInteraction
            String mitab = response.getResultSet().getMitab();

            // job exists and some data are associated with it
            if (mitab != null){
                PsimiTabReader reader = new PsimiTabReader();
                final List<BinaryInteraction> binaryInteractions = new ArrayList<BinaryInteraction>(maxResults);

				Iterator<BinaryInteraction> iterator;

				iterator = reader.iterate(mitab);
				while (iterator.hasNext()) {
					BinaryInteraction interaction = iterator.next();
					//We use the original query only to flip, because we ask the index all the results (query is *);
					flipIfNecessary(interaction, originalQuery);
					binaryInteractions.add(interaction);
				}

                result = new SearchResult<BinaryInteraction>( binaryInteractions,
                        response.getResultInfo().getTotalResults(),
                        firstResult,
                        maxResults,
                        null );
            }
            else {
                result = new SearchResult<BinaryInteraction>( Collections.EMPTY_LIST,
                        response.getResultInfo().getTotalResults(),
                        firstResult,
                        maxResults,
                        null );
            }
        } catch ( IOException e ) {
            throw new IllegalStateException( "Could not parse clustered MITAB data", e );
        } catch ( ClusteringServiceException e ) {
            throw new IllegalStateException( "Could not query clustered service", e );
        }
    }

	private void flipIfNecessary(BinaryInteraction interaction, String searchQuery) {
		flipInteractors.flipIfNecessary(interaction, searchQuery);
	}

	@Override
    public int getRowCount() {
        return Long.valueOf(result.getTotalCount()).intValue();
    }

	@Override
	public Object getRowData() {
        if (result == null) {
            return null;
        }

        if (!isRowWithinResultRange()) {
            firstResult = getRowIndex();
            fetchResults();
        }

        if (!isRowAvailable()) {
			FacesContext context = FacesContext.getCurrentInstance();
			String msg = "Temporarily impossible to retrieve results";
			String details = "Error while building results based on the query: '" + query + "'";

			if (context != null) {
				FacesMessage facesMessage = new FacesMessage(FacesMessage.SEVERITY_ERROR, msg, details);
				context.addMessage(null, facesMessage);
			}
			log.error("Row is unavailable");
			return null;
//            throw new IllegalArgumentException("row is unavailable");
        }

        final BinaryInteraction binaryInteraction = result.getData().get(rowIndex - firstResult);
        return binaryInteraction;
    }

	@Override
	public int getRowIndex() {
        return rowIndex;
    }

	@Override
	public Object getWrappedData() {
        return result;
    }

	@Override
	public boolean isRowAvailable() {
        if (result == null) {
            return false;
        }

        return rowIndex >= 0 && rowIndex < result.getTotalCount();
    }

    protected boolean isRowWithinResultRange() {
        return (getRowIndex() >= firstResult) && (getRowIndex() < (firstResult+maxResults));
    }

	@Override
	public void setRowIndex(int rowIndex) {
        if (rowIndex < -1) {
            throw new IllegalArgumentException("illegal rowIndex " + rowIndex);
        }
        int oldRowIndex = rowIndex;
        this.rowIndex = rowIndex;
        if (result != null && oldRowIndex != this.rowIndex) {
            Object data = isRowAvailable() ? getRowData() : null;
            DataModelEvent event = new DataModelEvent(this, this.rowIndex, data);
            DataModelListener[] listeners = getDataModelListeners();
            for (int i = 0; i < listeners.length; i++) {
                listeners[i].rowSelected(event);
            }
        }
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
