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
package org.hupo.psi.mi.psicquic.view.webapp.controller.search;

import org.apache.myfaces.orchestra.conversation.annotations.ConversationName;
import org.hupo.psi.mi.psicquic.view.webapp.controller.BaseController;
import org.hupo.psi.mi.psicquic.view.webapp.controller.config.PsicquicViewConfig;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import javax.annotation.PostConstruct;
import javax.faces.event.ActionEvent;
import javax.faces.model.SelectItem;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * User query object wrapper.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
@Controller("userQuery")
@Scope("conversation.access")
@ConversationName("general")
public class UserQuery extends BaseController {

    public static final String STAR_QUERY = "*:*";

    @Autowired
    private PsicquicViewConfig config;

    private String searchQuery = "*";

    private Map<String, String> termMap = new HashMap<String, String>();

    //for sorting and ordering
    private static final String DEFAULT_SORT_COLUMN = "rigid";
    private static final boolean DEFAULT_SORT_ORDER = true;

    private String userSortColumn = DEFAULT_SORT_COLUMN;
    private boolean userSortOrder = DEFAULT_SORT_ORDER;

    private int pageSize = 25;

    private boolean showNewFieldPanel;
    private QueryToken newQueryToken;


    private SearchField[] searchFields;

    private List<SelectItem> searchFieldSelectItems;
    private Map<String, SearchField> searchFieldsMap;

    public UserQuery() {
    }

    @PostConstruct
    public void reset() {
        this.searchQuery = "*";
        this.userSortColumn = DEFAULT_SORT_COLUMN;
        this.userSortOrder = DEFAULT_SORT_ORDER;

        clearFilters();

        initSearchFields();


    }

    public void clearFilters() {
        termMap.clear();
    }

    private void initSearchFields() {

        searchFields = new SearchField[]{
                new SearchField("", "All"),
                new SearchField("id", "Interactor id (Ex: P74565)"),
                new SearchField("alias", "Interactor alias (Ex: KHDRBS1)"),
                new SearchField("identifier", "Interactor id or alias"),
                new SearchField("pubauth", "Author (Ex: scott)"),
                new SearchField("pubid", "Publication id (Ex: 10837477)"),
                new SearchField("species", "Organism (Ex: human)"),
                new SearchField("type", "Interaction type (Ex: physical association)"),
                new SearchField("detmethod", "Interaction detection method (Ex: pull down)"),
                new SearchField("interaction_id", "Interaction Id (Ex: EBI-761050)"),
                new SearchField("pbiorole", "Biological roles (Ex: enzyme)"),
                new SearchField("ptype", "Interactor type (Ex: protein)"),
                new SearchField("pxref", "Interactor xref (Ex: GO:0003824)"),
                new SearchField("xref", "Interaction xref (Ex: nuclear pore)"),
                new SearchField("annot", "Interaction annotation (Ex: imex curation)"),
                new SearchField("udate", "Last update date (Ex: [YYYYMMDD TO YYYYMMDD])"),
                new SearchField("negative", "Negative interaction (Ex: true)", listNegativeSelectItems()),
                new SearchField("complex", "Complex expansion (Ex: spoke)"),
                new SearchField("ftype", "Interactor feature (Ex: binding site)"),
                new SearchField("pmethod", "Interactor identification method (Ex: western blot)"),
                new SearchField("stc", "Stoichiometry (Ex: true)", listStoichiometrySelectItems()),
                new SearchField("param", "Interaction parameter (Ex: true)", listParametersSelectItems())
		};


        searchFieldSelectItems = new ArrayList<SelectItem>(searchFields.length);
        searchFieldsMap = new HashMap<String, SearchField>();

        for (SearchField searchField : searchFields) {
            searchFieldSelectItems.add(new SelectItem(searchField.getName(), searchField.getDisplayName()));
            searchFieldsMap.put(searchField.getName(), searchField);
        }

   }

	private List<SelectItem> listNegativeSelectItems() {
		List<SelectItem> negativeItems = new ArrayList<SelectItem>(3);
		negativeItems.add(new SelectItem("true", "Only negative interactions"));
		negativeItems.add(new SelectItem("(true OR false)", "Include negative interactions"));
		negativeItems.add(new SelectItem("false", "Only positive interactions (default)"));
		return negativeItems;
	}

	private List<SelectItem> listParametersSelectItems() {
		List<SelectItem> parameters = new ArrayList<SelectItem>(2);
		parameters.add(new SelectItem("true", "Only interactions having parameters available"));
		parameters.add(new SelectItem("false", "Excludes interactions having parameters available"));
		return parameters;
	}

	private List<SelectItem> listStoichiometrySelectItems() {
		List<SelectItem> stoichiometry = new ArrayList<SelectItem>(2);
		stoichiometry.add(new SelectItem("true", "Only interactions having stoichiometry information"));
		stoichiometry.add(new SelectItem("false", "Excludes interactions having stoichiometry information"));
		return stoichiometry;
	}

	public void clearSearchFilters(ActionEvent evt) {
        clearFilters();
    }


    public void doShowAddFieldPanel(ActionEvent evt) {
        showAddFieldsPanel();
        newQueryToken = new QueryToken("");
    }

    public void doAddFieldToQuery(ActionEvent evt) {
        doAddFieldToQuery(newQueryToken);
    }

    public void doAddFieldToQuery(QueryToken queryToken) {
        if (!isWildcardQuery(queryToken.getQuery())) {

            if (isWildcardQuery(searchQuery)) {
                searchQuery = queryToken.toQuerySyntax(true);
            } else {
                searchQuery = surroundByBracesIfNecessary(searchQuery);
                searchQuery = searchQuery + " " + queryToken.toQuerySyntax();
            }
        }

        hideAddFieldsPanel();
    }

    public boolean isWildcardQuery() {
        return isWildcardQuery(searchQuery);
    }

    private boolean isWildcardQuery(String query) {
        return (query == null || query.trim().length() == 0 || "*".equals(query) || STAR_QUERY.equals(query));
    }

    private String surroundByBracesIfNecessary(String query) {
        if (query.contains(" AND ") || query.contains(" OR ")) {
            query = "(" + query + ")";
        }

        return query;
    }

    public void doCancelAddField(ActionEvent evt) {
        hideAddFieldsPanel();
    }

    private void showAddFieldsPanel() {
        showNewFieldPanel = true;
    }

    private void hideAddFieldsPanel() {
        newQueryToken = null;
        showNewFieldPanel = false;
    }

    public String getSearchQuery() {
        if ("".equals(searchQuery) || STAR_QUERY.equals(searchQuery)) {
            searchQuery = "*";
        }
        return searchQuery;
    }

    public String getFilteredSearchQuery() {
        String query = getSearchQuery();
        String miqlFilterQuery = config.getMiqlFilterQuery();

        if (miqlFilterQuery == null || miqlFilterQuery.length() == 0) {
            return query;
        } else {
            if (query.equals("*") || query.length() == 0) {
                query = miqlFilterQuery;
            } else {
                query = "(" + query + ") AND (" + miqlFilterQuery + ")";
            }
        }

        return query;

    }

    public void setSearchQuery(String searchQuery) {
        this.searchQuery = searchQuery;
    }

    public void resetSearchQuery() {
        this.searchQuery = null;
    }

    public String getUserSortColumn() {
        return userSortColumn;
    }

    public void setUserSortColumn(String userSortColumn) {
        this.userSortColumn = userSortColumn;
    }

    public boolean getUserSortOrder() {
        return userSortOrder;
    }

    public void setUserSortOrder(boolean userSortOrder) {
        this.userSortOrder = userSortOrder;
    }

    public int getPageSize() {
        return pageSize;
    }

    public void setPageSize(int pageSize) {
        this.pageSize = pageSize;
    }

    public Map<String, String> getTermMap() {
        return termMap;
    }

    public void setTermMap(Map<String, String> termMap) {
        this.termMap = termMap;
    }

    public List<SelectItem> getSearchFieldSelectItems() {
        return searchFieldSelectItems;
    }

    public boolean isShowNewFieldPanel() {
        return showNewFieldPanel;
    }

    public QueryToken getNewQueryToken() {
        return newQueryToken;
    }

    public void setNewQueryToken(QueryToken newQueryToken) {
        this.newQueryToken = newQueryToken;
    }

    public Map<String, SearchField> getSearchFieldsMap() {
        return searchFieldsMap;
    }
}
