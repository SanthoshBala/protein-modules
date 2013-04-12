package org.hupo.psi.mi.psicquic.model;

import org.apache.solr.client.solrj.response.FacetField;
import org.apache.solr.common.SolrDocumentList;
import org.hupo.psi.calimocho.io.IllegalRowException;
import org.hupo.psi.calimocho.tab.model.ColumnBasedDocumentDefinition;
import org.hupo.psi.calimocho.tab.util.MitabDocumentDefinitionFactory;
import org.hupo.psi.calimocho.xgmml.XGMMLDocumentDefinition;
import org.hupo.psi.mi.rdf.PsimiRdfConverter;
import org.hupo.psi.mi.rdf.RdfFormat;
import psidev.psi.mi.tab.PsimiTabException;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.converter.tab2xml.Tab2Xml;
import psidev.psi.mi.tab.converter.tab2xml.XmlConversionException;
import psidev.psi.mi.xml.PsimiXmlVersion;
import psidev.psi.mi.xml.model.EntrySet;

import java.io.*;
import java.util.Collections;
import java.util.List;

/**
 * Contains the results of a psicquic search. It can create XML, XGMML or any RDF/Biopax format from these results
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>27/07/12</pre>
 */

public class PsicquicSearchResults {

    private long numberResults;
    private SolrDocumentList results;
    private List<FacetField> facetFields;
    private String [] fieldNames;

    private PsimiTabReader mitabReader;
    private PsimiRdfConverter rdfConverter;

    public PsicquicSearchResults(SolrDocumentList results, String [] fieldNames){
        this.results = results;
        this.fieldNames = fieldNames;
        mitabReader = new PsimiTabReader();
        rdfConverter = new PsimiRdfConverter();

        if (results != null){
            this.numberResults = results.getNumFound();
        }
        else {
            this.numberResults = 0;
        }
    }

    public PsicquicSearchResults(SolrDocumentList results, String [] fieldNames, List<FacetField> facetFields){
        this.results = results;
        this.fieldNames = fieldNames;
        mitabReader = new PsimiTabReader();
        rdfConverter = new PsimiRdfConverter();

        if (results != null){
            this.numberResults = results.getNumFound();
        }
        else {
            this.numberResults = 0;
        }

        this.facetFields = facetFields;
    }

    public InputStream getMitab() {
        if (this.fieldNames == null || this.fieldNames.length == 0){
            return null;
        }

        return new MitabInputStream(this.results, this.fieldNames);
    }

    public long getNumberResults() {
        return numberResults;
    }

    public EntrySet createEntrySet() throws PsimiTabException, IOException, XmlConversionException, IllegalAccessException {
        if (this.numberResults == 0 || this.results.size() == 0) {
            return new EntrySet();
        }

        Tab2Xml tab2Xml = new Tab2Xml();
        InputStream inputStream = getMitab();
        EntrySet entrySet = null;
        try{
            if (inputStream != null){
                entrySet = tab2Xml.convert(mitabReader.read(inputStream));
            }
            else {
                entrySet = new EntrySet(PsimiXmlVersion.VERSION_25_UNDEFINED);
            }
        }
        finally {
            if (inputStream != null){
                inputStream.close();
            }
            tab2Xml.close();
        }
        return entrySet;
    }

    public InputStream createXGMML() throws IOException, IllegalRowException {
        InputStream mitabOs = getMitab();
        InputStream xgmml = null;

        try{
            ColumnBasedDocumentDefinition mitab25Definition = MitabDocumentDefinitionFactory.mitab25();
            XGMMLDocumentDefinition xgmmlDefinition = new XGMMLDocumentDefinition("PSICQUIC", "Generated from MITAB 2.5", "http://psicquic.googlecode.com");

            ByteArrayOutputStream xgmmlOutput = new ByteArrayOutputStream();
            Writer xgmmlWriter = new BufferedWriter(new OutputStreamWriter(xgmmlOutput));

            try{

                xgmmlDefinition.writeDocument(xgmmlWriter, mitabOs, mitab25Definition);
                xgmml = new ByteArrayInputStream(xgmmlOutput.toByteArray());

            }
            finally {
                // close stringWriter now
                xgmmlWriter.close();
                xgmmlOutput.close();
            }
        }
        finally {
            // close mitabOs now
            if (mitabOs != null){
                mitabOs.close();
            }
        }

        return xgmml;
    }

    protected String getRdfFormatName(String format) {
        if (format.equalsIgnoreCase("biopax") || format.equalsIgnoreCase("biopax-L3")) {
            return RdfFormat.BIOPAX_L3.getName();
        } else if (format.equalsIgnoreCase("biopax-L2")) {
            return RdfFormat.BIOPAX_L2.getName();
        }

        format = format.substring(4);

        String rdfFormat;

        if ("xml".equalsIgnoreCase(format)) {
            rdfFormat = "RDF/XML";
        } else if ("xml-abbrev".equalsIgnoreCase(format)) {
            rdfFormat = "RDF/XML-ABBREV";
        } else {
            rdfFormat = format.toUpperCase();
        }

        return rdfFormat;
    }

    public InputStream createRDFOrBiopax(String format) throws PsimiTabException, IOException, XmlConversionException, IllegalAccessException {

        if (format == null || format.length() < 4){
            return null;
        }

        String rdfFormat = getRdfFormatName(format);

        psidev.psi.mi.xml.model.EntrySet entrySet = createEntrySet();

        ByteArrayOutputStream rdfOutput = new ByteArrayOutputStream();
        Writer rdfWriter = new BufferedWriter(new OutputStreamWriter(rdfOutput));
        InputStream rdf = null;

        try{
            // only convert when having some results
            if (entrySet != null && !entrySet.getEntries().isEmpty()){
                rdfConverter.convert(entrySet, rdfFormat , rdfWriter);
                rdf = new ByteArrayInputStream(rdfOutput.toByteArray());
            }
        }
        finally {
            // close writer
            rdfWriter.close();
            rdfOutput.close();
            rdfConverter.close();
        }

        return rdf;
    }

    /**
     *
     * @return the solr document list
     */
    public SolrDocumentList getSolrDocumentList() {
        return results;
    }

    /**
     *
     * @return the list of facet fields that were matching query
     */
    public List<FacetField> getFacetFieldList() {
        return facetFields != null ? facetFields : Collections.EMPTY_LIST;
    }

    /**
     *
     * @param facetFieldName
     * @return the facet field with a specific name, null if it does not exist
     */
    public FacetField getFacetField(String facetFieldName){
        if (this.facetFields==null || facetFieldName == null) return null;
        for (FacetField f : this.facetFields) {
            if (f.getName().equals(facetFieldName)) return f;
        }
        return null;
    }
}
