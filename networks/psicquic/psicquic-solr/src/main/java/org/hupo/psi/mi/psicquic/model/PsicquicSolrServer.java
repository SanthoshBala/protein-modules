package org.hupo.psi.mi.psicquic.model;

import org.apache.commons.lang.ArrayUtils;
import org.apache.solr.client.solrj.SolrQuery;
import org.apache.solr.client.solrj.SolrServer;
import org.apache.solr.client.solrj.SolrServerException;
import org.apache.solr.client.solrj.impl.HttpSolrServer;
import org.apache.solr.client.solrj.response.FacetField;
import org.apache.solr.common.SolrDocumentList;
import org.apache.solr.common.params.FacetParams;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import psidev.psi.mi.calimocho.solr.converter.SolrFieldName;
import psidev.psi.mi.tab.PsimiTabReader;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Psicquic solr server that is wrapping a solrServer
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>12/06/12</pre>
 */

public class PsicquicSolrServer {

    private final Logger logger = LoggerFactory.getLogger(PsicquicSolrServer.class);

    // static final variables
    protected final static String STORED_FIELD_EXTENSION="_o";
    public static final String RETURN_TYPE_XML25 = "psi-mi/xml25";
    public static final String RETURN_TYPE_MITAB25 = "psi-mi/tab25";
    public static final String RETURN_TYPE_MITAB26 = "psi-mi/tab26";
    public static final String RETURN_TYPE_MITAB27 = "psi-mi/tab27";
    public static final String RETURN_TYPE_COUNT = "count";
    public static final String RETURN_TYPE_DEFAULT = RETURN_TYPE_MITAB25;
    private final static String DISMAX_PARAM_NAME = "qf";
    private final static String DEFAULT_PARAM_NAME = "df";
    private final static String DEFAULT_MM_PARAM = "mm";
    private final static String DISMAX_TYPE = "edismax";
    private final static String QUERY_TYPE = "defType";

    /**
     * solr server
     */
    protected SolrServer solrServer;
    /**
     * MITAB reader
     */
    protected PsimiTabReader mitabReader;
    /**
     * Map of SOLR fields that can be returned for a specific format
     */
    protected Map<String, String []> solrFields;

    public static String[] DATA_FIELDS_25 = new String[] {
            SolrFieldName.idA+STORED_FIELD_EXTENSION, SolrFieldName.idB+STORED_FIELD_EXTENSION, SolrFieldName.altidA+STORED_FIELD_EXTENSION,
            SolrFieldName.altidB+STORED_FIELD_EXTENSION, SolrFieldName.aliasA+STORED_FIELD_EXTENSION, SolrFieldName.aliasB+STORED_FIELD_EXTENSION,
            SolrFieldName.detmethod+STORED_FIELD_EXTENSION, SolrFieldName.pubauth+STORED_FIELD_EXTENSION, SolrFieldName.pubid+STORED_FIELD_EXTENSION,
            SolrFieldName.taxidA+STORED_FIELD_EXTENSION, SolrFieldName.taxidB+STORED_FIELD_EXTENSION, SolrFieldName.type+STORED_FIELD_EXTENSION,
            SolrFieldName.source+STORED_FIELD_EXTENSION,SolrFieldName.interaction_id+STORED_FIELD_EXTENSION, SolrFieldName.confidence+STORED_FIELD_EXTENSION
    };

    public static String[] DATA_FIELDS_26 = new String[] {
            SolrFieldName.idA+STORED_FIELD_EXTENSION, SolrFieldName.idB+STORED_FIELD_EXTENSION, SolrFieldName.altidA+STORED_FIELD_EXTENSION,
            SolrFieldName.altidB+STORED_FIELD_EXTENSION, SolrFieldName.aliasA+STORED_FIELD_EXTENSION, SolrFieldName.aliasB+STORED_FIELD_EXTENSION,
            SolrFieldName.detmethod+STORED_FIELD_EXTENSION, SolrFieldName.pubauth+STORED_FIELD_EXTENSION, SolrFieldName.pubid+STORED_FIELD_EXTENSION,
            SolrFieldName.taxidA+STORED_FIELD_EXTENSION, SolrFieldName.taxidB+STORED_FIELD_EXTENSION, SolrFieldName.type+STORED_FIELD_EXTENSION,
            SolrFieldName.source+STORED_FIELD_EXTENSION,SolrFieldName.interaction_id+STORED_FIELD_EXTENSION, SolrFieldName.confidence+STORED_FIELD_EXTENSION,
            SolrFieldName.complex+STORED_FIELD_EXTENSION, SolrFieldName.pbioroleA+STORED_FIELD_EXTENSION, SolrFieldName.pbioroleB+STORED_FIELD_EXTENSION,
            SolrFieldName.pexproleA+STORED_FIELD_EXTENSION, SolrFieldName.pexproleB+STORED_FIELD_EXTENSION, SolrFieldName.ptypeA+STORED_FIELD_EXTENSION,
            SolrFieldName.ptypeB+STORED_FIELD_EXTENSION, SolrFieldName.pxrefA+STORED_FIELD_EXTENSION, SolrFieldName.pxrefB+STORED_FIELD_EXTENSION,
            SolrFieldName.xref+STORED_FIELD_EXTENSION, SolrFieldName.annotA+STORED_FIELD_EXTENSION,
            SolrFieldName.annotB+STORED_FIELD_EXTENSION, SolrFieldName.annot+STORED_FIELD_EXTENSION, SolrFieldName.taxidHost+STORED_FIELD_EXTENSION,
            SolrFieldName.param+STORED_FIELD_EXTENSION, SolrFieldName.cdate+STORED_FIELD_EXTENSION, SolrFieldName.udate+STORED_FIELD_EXTENSION,
            SolrFieldName.checksumA+STORED_FIELD_EXTENSION, SolrFieldName.checksumB+STORED_FIELD_EXTENSION, SolrFieldName.checksumI+STORED_FIELD_EXTENSION,
            SolrFieldName.negative+STORED_FIELD_EXTENSION
    };

    public static String[] DATA_FIELDS_27 = new String[] {
            SolrFieldName.idA+STORED_FIELD_EXTENSION, SolrFieldName.idB+STORED_FIELD_EXTENSION, SolrFieldName.altidA+STORED_FIELD_EXTENSION,
            SolrFieldName.altidB+STORED_FIELD_EXTENSION, SolrFieldName.aliasA+STORED_FIELD_EXTENSION, SolrFieldName.aliasB+STORED_FIELD_EXTENSION,
            SolrFieldName.detmethod+STORED_FIELD_EXTENSION, SolrFieldName.pubauth+STORED_FIELD_EXTENSION, SolrFieldName.pubid+STORED_FIELD_EXTENSION,
            SolrFieldName.taxidA+STORED_FIELD_EXTENSION, SolrFieldName.taxidB+STORED_FIELD_EXTENSION, SolrFieldName.type+STORED_FIELD_EXTENSION,
            SolrFieldName.source+STORED_FIELD_EXTENSION,SolrFieldName.interaction_id+STORED_FIELD_EXTENSION, SolrFieldName.confidence+STORED_FIELD_EXTENSION,
            SolrFieldName.complex+STORED_FIELD_EXTENSION, SolrFieldName.pbioroleA+STORED_FIELD_EXTENSION, SolrFieldName.pbioroleB+STORED_FIELD_EXTENSION,
            SolrFieldName.pexproleA+STORED_FIELD_EXTENSION, SolrFieldName.pexproleB+STORED_FIELD_EXTENSION, SolrFieldName.ptypeA+STORED_FIELD_EXTENSION,
            SolrFieldName.ptypeB+STORED_FIELD_EXTENSION, SolrFieldName.pxrefA+STORED_FIELD_EXTENSION, SolrFieldName.pxrefB+STORED_FIELD_EXTENSION,
            SolrFieldName.xref+STORED_FIELD_EXTENSION, SolrFieldName.annotA+STORED_FIELD_EXTENSION,
            SolrFieldName.annotB+STORED_FIELD_EXTENSION, SolrFieldName.annot+STORED_FIELD_EXTENSION, SolrFieldName.taxidHost+STORED_FIELD_EXTENSION,
            SolrFieldName.param+STORED_FIELD_EXTENSION, SolrFieldName.cdate+STORED_FIELD_EXTENSION, SolrFieldName.udate+STORED_FIELD_EXTENSION,
            SolrFieldName.checksumA+STORED_FIELD_EXTENSION, SolrFieldName.checksumB+STORED_FIELD_EXTENSION, SolrFieldName.checksumI+STORED_FIELD_EXTENSION,
            SolrFieldName.negative+STORED_FIELD_EXTENSION, SolrFieldName.ftypeA+STORED_FIELD_EXTENSION, SolrFieldName.ftypeB+STORED_FIELD_EXTENSION,
            SolrFieldName.stcA+STORED_FIELD_EXTENSION, SolrFieldName.stcB+STORED_FIELD_EXTENSION, SolrFieldName.pmethodA+STORED_FIELD_EXTENSION,
            SolrFieldName.pmethodB+STORED_FIELD_EXTENSION
    };

    public PsicquicSolrServer(SolrServer solrServer){
        this.solrServer = solrServer;

        if (this.solrServer == null){
            throw new IllegalArgumentException("Cannot create a new PsicquicSolrServer if the SolrServer is null");
        }

        // initialise default solr field map
        initializeSolrFieldsMap();

        // create new mitab reader without header
        mitabReader = new PsimiTabReader();
    }

    /**
     * Associates for each format the fields that are expected to be returned
     */
    protected void initializeSolrFieldsMap(){
        solrFields = new HashMap<String, String[]>(5);
        solrFields.put(RETURN_TYPE_MITAB25, DATA_FIELDS_25);
        solrFields.put(RETURN_TYPE_MITAB26, DATA_FIELDS_26);
        solrFields.put(RETURN_TYPE_MITAB27, DATA_FIELDS_27);
        solrFields.put(RETURN_TYPE_XML25, DATA_FIELDS_27);
        solrFields.put(RETURN_TYPE_COUNT, new String[] {});
    }

    /**
     *
     * @param q : query entered by user
     * @param firstResult : first result
     * @param maxResults : max number of results
     * @param returnType : format of the results
     * @param queryFilter : any filter to add to the query
     * @return  the response
     * @throws PsicquicSolrException
     */
    public PsicquicSearchResults search(String q, Integer firstResult, Integer maxResults, String returnType, String queryFilter) throws PsicquicSolrException, SolrServerException {

        return searchWithFilters(q, firstResult, maxResults, returnType, queryFilter != null ? new String[]{queryFilter} : null);
    }

    public PsicquicSearchResults searchWithFilters(String q, Integer firstResult, Integer maxResults, String returnType, String [] queryFilter) throws PsicquicSolrException, SolrServerException {
        if (q == null) throw new NullPointerException("Null query");

        // format wildcard query
        if ("*".equals(q)){
            q = "*:*";
        }

        SolrQuery solrQuery = new SolrQuery(q);

        // use dismax parser for querying default fields
        //solrQuery.setParam(DISMAX_PARAM_NAME, SolrFieldName.identifier.toString(), SolrFieldName.pubid.toString(), SolrFieldName.pubauth.toString(), SolrFieldName.species.toString(), SolrFieldName.detmethod.toString(), SolrFieldName.type.toString(), SolrFieldName.interaction_id.toString());
        solrQuery.setParam(DISMAX_PARAM_NAME, SolrFieldName.identifier.toString()+" "+SolrFieldName.pubid.toString()+" "
                +SolrFieldName.pubauth.toString()+" "+SolrFieldName.species.toString()+" "+SolrFieldName.detmethod.toString()
                +" "+SolrFieldName.type.toString()+" "+SolrFieldName.interaction_id.toString()+" "+SolrFieldName.xref.toString()
                +" "+SolrFieldName.pxref.toString());
        solrQuery.setParam(QUERY_TYPE, DISMAX_TYPE);
        solrQuery.setParam(DEFAULT_MM_PARAM, "1");

        // set first result
        if (firstResult != null)
        {
            solrQuery.setStart(firstResult);
        }else {
            solrQuery.setStart(0);
        }

        // set max results
        // WARNING in solr 3.6
        // * an *NumberFormatException* occurs if _rows_ > 2147483647
        // * an *ArrayIndexOutOfBoundsException* occurs if _rows_ + _start_ > 2147483647; e.g. _rows_ = 2147483640 and _start_ = 8
        // we need to substract to avoid this exception
        if (maxResults != null) {
            solrQuery.setRows(maxResults);
        } else {
            solrQuery.setRows(Integer.MAX_VALUE - solrQuery.getStart());
        }

        // by default, if no return types is specified, it will return MITAB 2.5
        if (returnType == null){
            returnType = RETURN_TYPE_DEFAULT;
        }

        // apply any filter
        if (queryFilter != null && queryFilter.length > 0) {
            for (String filter : queryFilter){
                if (!"*".equals(filter)) {
                    solrQuery.addFilterQuery(filter);
                }
            }
        }

        // use normal query so we can use dismax parser to define multiple default field names in case of free text searches. The filter query fq parameter does not work with
        // dismax parser
        solrQuery.setQuery(q);

        return search(solrQuery, returnType);
    }

    /**
     * Allows faceting search
     * @param q : main user query
     * @param firstResult : first result of the query
     * @param maxResults : maximun number of results returned by the query
     * @param returnType : type of the results
     * @param queryFilter : query filters
     * @param facets : names of the facet fields to return with this query
     * @param firstFacet : first facet elements in the facet field lists (for pagination)
     * @param maxFacet : maxNumber of facet fields to return
     * @return
     * @throws PsicquicSolrException
     * @throws SolrServerException
     */
    public PsicquicSearchResults searchWithFacets(String q, Integer firstResult, Integer maxResults, String returnType, String [] queryFilter, String [] facets, Integer firstFacet, Integer maxFacet) throws PsicquicSolrException, SolrServerException {
        if (q == null) throw new NullPointerException("Null query");

        // format wildcard query
        if ("*".equals(q)){
            q = "*:*";
        }

        SolrQuery solrQuery = new SolrQuery(q);

        // use dismax parser for querying default fields
        //solrQuery.setParam(DISMAX_PARAM_NAME, SolrFieldName.identifier.toString(), SolrFieldName.pubid.toString(), SolrFieldName.pubauth.toString(), SolrFieldName.species.toString(), SolrFieldName.detmethod.toString(), SolrFieldName.type.toString(), SolrFieldName.interaction_id.toString());
        solrQuery.setParam(DISMAX_PARAM_NAME, SolrFieldName.identifier.toString()+" "+SolrFieldName.pubid.toString()+" "
                +SolrFieldName.pubauth.toString()+" "+SolrFieldName.species.toString()+" "+SolrFieldName.detmethod.toString()
                +" "+SolrFieldName.type.toString()+" "+SolrFieldName.interaction_id.toString()+" "+SolrFieldName.xref.toString()
                +" "+SolrFieldName.pxref.toString());        solrQuery.setParam(QUERY_TYPE, DISMAX_TYPE);
        solrQuery.setParam(DEFAULT_MM_PARAM, "1");

        // set first result
        if (firstResult != null)
        {
            solrQuery.setStart(firstResult);
        }else {
            solrQuery.setStart(0);
        }

        // set max results
        // WARNING in solr 3.6
        // * an *NumberFormatException* occurs if _rows_ > 2147483647
        // * an *ArrayIndexOutOfBoundsException* occurs if _rows_ + _start_ > 2147483647; e.g. _rows_ = 2147483640 and _start_ = 8
        // we need to substract to avoid this exception
        if (maxResults != null) {
            solrQuery.setRows(maxResults);
        } else {
            solrQuery.setRows(Integer.MAX_VALUE - solrQuery.getStart());
        }

        // allows faceting
        if (facets != null && facets.length > 0){
            // we allow faceting
            solrQuery.setFacet(true);

            // we want all the facet fields with min count = 1. The facet fields with count = 0 are not interesting
            solrQuery.setFacetMinCount(1);

            // set first result
            if (firstFacet != null)
            {
                solrQuery.set(FacetParams.FACET_OFFSET, firstFacet);
            }else {
                firstFacet = 0;
                solrQuery.set(FacetParams.FACET_OFFSET, 0);
            }

            // set max results
            // WARNING in solr 3.6
            // * an *NumberFormatException* occurs if _rows_ > 2147483647
            // * an *ArrayIndexOutOfBoundsException* occurs if _rows_ + _start_ > 2147483647; e.g. _rows_ = 2147483640 and _start_ = 8
            // we need to substract to avoid this exception
            if (maxFacet != null) {
                solrQuery.setFacetLimit(maxFacet);
            } else {
                solrQuery.setFacetLimit(Integer.MAX_VALUE - firstFacet);
            }

            solrQuery.setFacetSort(FacetParams.FACET_SORT_COUNT);

            // add facet fields
            solrQuery.addFacetField(facets);

        }

        // by default, if no return types is specified, it will return MITAB 2.5
        if (returnType == null){
            returnType = RETURN_TYPE_DEFAULT;
        }

        // apply any filter
        if (queryFilter != null && queryFilter.length > 0) {
            for (String filter : queryFilter){
                if (!"*".equals(filter)) {
                    solrQuery.addFilterQuery(filter);
                }
            }
        }

        // use normal query so we can use dismax parser to define multiple default field names in case of free text searches. The filter query fq parameter does not work with
        // dismax parser
        solrQuery.setQuery(q);

        return search(solrQuery, returnType);
    }

    protected PsicquicSearchResults search(SolrQuery originalQuery, String returnType) throws PsicquicSolrException, SolrServerException {

        String[] fields = solrFields.containsKey(returnType) ? solrFields.get(returnType) : new String[]{};
        SolrQuery copy = originalQuery.getCopy();

        if(copy.getFields()!=null){
            // new array created here so the map solrFields is never changed
            fields = (String[]) ArrayUtils.addAll(fields, copy.getFields().split(","));
        }

        copy.setFields(fields);

        String query = copy.getQuery();
        // if using a wildcard query we convert to lower case
        // as of http://mail-archives.apache.org/mod_mbox/lucene-solr-user/200903.mbox/%3CFD3AFB65-AEC1-40B2-A0A4-7E14A519AB05@ehatchersolutions.com%3E
        if (query.contains("*")) {
            String[] tokens = query.split(" ");

            StringBuilder sb = new StringBuilder(query.length());

            for (String token : tokens) {
                if (token.contains("*")) {
                    // name of the field with upper case to keep
                    if (token.startsWith("idA:") || token.startsWith("idB:") || token.startsWith("taxidA:") || token.startsWith("taxidB:")
                            || token.startsWith("pbioroleA:") || token.startsWith("pbioroleB:") || token.startsWith("ptypeA:")
                            || token.startsWith("ptypeB:") || token.startsWith("pxrefA:") || token.startsWith("pxrefB:")
                            || token.startsWith("ftypeA:") || token.startsWith("ftypeB:") || token.startsWith("pmethodA:") || token.startsWith("pmethodA:")){

                        int firstIndex = token.indexOf(":");
                        String prefix = token.substring(0, firstIndex+1);
                        String correctedValue = token.substring(firstIndex+1).toLowerCase();
                        sb.append(prefix).append(correctedValue);
                    }
                    else {
                        sb.append(token.toLowerCase());
                    }
                } else {
                    sb.append(token);
                }

                sb.append(" ");
            }

            query = sb.toString().trim();
        }

        // replace and/or in lower case
        copy.setQuery(query
                .replaceAll(" and ", " AND ")
                .replaceAll(" or ", " OR ")
                .replaceAll(" not ", " NOT ")
                .replaceAll("^not ", "NOT "));

        // By default filter negative interactions. If a neagtive field is specified in the query, then the filter which excludes negative
        // interactions is removed
        if (!containsNegativeFilter(copy) || returnType == null || (returnType != null && RETURN_TYPE_MITAB25.equals(returnType))){
            copy.addFilterQuery(SolrFieldName.negative+":false");
        }

        org.apache.solr.client.solrj.response.QueryResponse solrResponse = solrServer.query(copy);

        if (solrResponse == null){
            return null;
        }

        return createSearchResults(solrResponse.getResults(), returnType, solrResponse.getFacetFields());
    }

    private boolean containsNegativeFilter(SolrQuery query){
        if (query.getQuery() != null && query.getQuery().contains(SolrFieldName.negative+":")){
            return true;
        }

        String[] queries = query.getFilterQueries();

        if (queries == null || queries.length == 0){
            return false;
        }

        for (String filter : queries){
            if (filter.contains(SolrFieldName.negative+":")){
                return true;
            }
        }

        return false;
    }

    /**
     *
     * @param docList
     * @param returnType : build the results based on return type. If it is not provided, MITAB 2.5 is the default return type
     * @param facetFields : list of facet fields returned by query response
     * @return PsicquicSearchResults built on MITAB
     * @throws PsicquicSolrException
     */
    protected PsicquicSearchResults createSearchResults(SolrDocumentList docList, String returnType, List<FacetField> facetFields) throws PsicquicSolrException {

        String resultType = returnType != null ? returnType : RETURN_TYPE_DEFAULT;
        PsicquicSearchResults results = createMitabResultsForType(docList, resultType, facetFields);

        return results;
    }

    /**
     *
     * @param docList : the list of solr documents returned by the query response
     * @param mitabType : build the results based on mitab type
     * @param facetFields : list of facet fields returned by query response
     * @return PsicquicSearchResults built on MITAB
     * @throws PsicquicSolrException
     */
    protected PsicquicSearchResults createMitabResultsForType(SolrDocumentList docList, String mitabType, List<FacetField> facetFields) throws PsicquicSolrException {

        String [] fieldNames = solrFields.get(mitabType);

        if (fieldNames == null){
            throw new PsicquicSolrException("The format " + mitabType + " is not a recognised MITAB format");
        }

        PsicquicSearchResults searchResults = new PsicquicSearchResults(docList, fieldNames, facetFields);
        return searchResults;
    }

    /**
     * Shutdown solr servers
     */
    public void shutdown(){

        if (this.solrServer != null && this.solrServer instanceof HttpSolrServer){
            HttpSolrServer httpsolrServer = (HttpSolrServer) solrServer;
            httpsolrServer.shutdown();
        }
    }
}
