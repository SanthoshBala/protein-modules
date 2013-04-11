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
package org.hupo.psi.mi.psicquic.ws;

import org.apache.commons.lang.StringUtils;
import org.hupo.psi.mi.psicquic.*;
import org.hupo.psi.mi.psicquic.ws.config.PsicquicConfig;
import org.hupo.psi.mi.psicquic.ws.utils.CompressedStreamingOutput;
import org.hupo.psi.mi.psicquic.ws.utils.PsicquicStreamingOutput;
import org.hupo.psi.mi.psicquic.ws.utils.XgmmlStreamingOutput;
import org.hupo.psi.mi.rdf.PsimiRdfConverter;
import org.hupo.psi.mi.rdf.RdfFormat;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import psidev.psi.mi.xml.converter.ConverterException;
import psidev.psi.mi.xml.converter.impl254.EntrySetConverter;
import psidev.psi.mi.xml.dao.inMemory.InMemoryDAOFactory;
import psidev.psi.mi.xml.io.impl.PsimiXmlWriter254;
import psidev.psi.mi.xml254.jaxb.Entry;
import psidev.psi.mi.xml254.jaxb.EntrySet;

import javax.ws.rs.core.GenericEntity;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.Map;


/**
 * This web service is based on a PSIMITAB SOLR index to search and return the results.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: IntactPsicquicService.java 12873 2009-03-18 02:51:31Z baranda $
 */
@Controller
public class IndexBasedPsicquicRestService implements PsicquicRestService {

    public static final String RETURN_TYPE_XML25 = "xml25";
    public static final String RETURN_TYPE_MITAB25 = "tab25";
    public static final String RETURN_TYPE_MITAB25_BIN = "tab25-bin";
    public static final String RETURN_TYPE_BIOPAX = "biopax";
    public static final String RETURN_TYPE_BIOPAX_l2 = "biopax-L2";
    public static final String RETURN_TYPE_BIOPAX_l3 = "biopax-L3";
    public static final String RETURN_TYPE_XGMML = "xgmml";
    public static final String RETURN_TYPE_RDF_XML = "rdf-xml";
    public static final String RETURN_TYPE_RDF_XML_ABBREV = "rdf-xml-abbrev";
    public static final String RETURN_TYPE_RDF_N3 = "rdf-n3";
    public static final String RETURN_TYPE_RDF_TURTLE = "rdf-turtle";
    public static final String RETURN_TYPE_COUNT = "count";

    @Autowired
    private PsicquicConfig config;

    @Autowired
    private PsicquicService psicquicService;

    public static final List<String> SUPPORTED_REST_RETURN_TYPES = Arrays.asList(
            RETURN_TYPE_XML25,
            RETURN_TYPE_MITAB25,
            RETURN_TYPE_BIOPAX,
            RETURN_TYPE_BIOPAX_l2,
            RETURN_TYPE_BIOPAX_l3,
            RETURN_TYPE_XGMML,
            RETURN_TYPE_RDF_XML,
            RETURN_TYPE_RDF_XML_ABBREV,
            RETURN_TYPE_RDF_N3,
            RETURN_TYPE_RDF_TURTLE,
            RETURN_TYPE_COUNT);

    public IndexBasedPsicquicRestService() {
    }

    public Response getByInteractor(String interactorAc, String db, String format, String firstResult, String maxResults, String compressed) throws PsicquicServiceException, NotSupportedMethodException, NotSupportedTypeException {
        String query = "identifier:"+createQueryValue(interactorAc, db);
        return getByQuery(query, format, firstResult, maxResults, compressed);
    }

    public Response getByInteraction(String interactionAc, String db, String format, String firstResult, String maxResults, String compressed) throws PsicquicServiceException, NotSupportedMethodException, NotSupportedTypeException {
        String query = "interaction_id:"+createQueryValue(interactionAc, db);
        return getByQuery(query, format, firstResult, maxResults, compressed);
    }

    public Response getByQuery(String query, String format,
                                                 String firstResultStr,
                                                 String maxResultsStr,
                                                 String compressed) throws PsicquicServiceException,
                                                                 NotSupportedMethodException,
                                                                 NotSupportedTypeException {

        boolean isCompressed = ("y".equalsIgnoreCase(compressed) || "true".equalsIgnoreCase(compressed));

        int firstResult;
        int maxResults;

        try {
            firstResult =Integer.parseInt(firstResultStr);
        } catch (NumberFormatException e) {
            throw new PsicquicServiceException("firstResult parameter is not a number: "+firstResultStr);
        }

        try {
            if (maxResultsStr == null) {
                maxResults = Integer.MAX_VALUE;
            } else {
                maxResults = Integer.parseInt(maxResultsStr);
            }
        } catch (NumberFormatException e) {
            throw new PsicquicServiceException("maxResults parameter is not a number: "+maxResultsStr);
        }

        format = format.toLowerCase();

        // if using mitab25-bin, set to mitab and compressed=y
        if (RETURN_TYPE_MITAB25_BIN.equalsIgnoreCase(format)) {
            format = RETURN_TYPE_MITAB25;
            isCompressed = true;
        }

        try {
            if (RETURN_TYPE_XML25.equalsIgnoreCase(format)) {
                final EntrySet entrySet = getByQueryXml(query, firstResult, maxResults);

                int count = 0;
                if (entrySet.getEntries() != null && !entrySet.getEntries().isEmpty()){
                    for (Entry entry : entrySet.getEntries()){
                        if (entry.getInteractionList() != null && entry.getInteractionList().getInteractions() != null && !entry.getInteractionList().getInteractions().isEmpty()){
                            count += entry.getInteractionList().getInteractions().size();
                        }
                    }
                }

                return prepareResponse(Response.status(200).type(MediaType.APPLICATION_XML), entrySet, count, isCompressed).build();
            } else if ((format.toLowerCase().startsWith("rdf") && format.length() > 5) || format.toLowerCase().startsWith("biopax")
                    || format.toLowerCase().startsWith("biopax-L3") || format.toLowerCase().startsWith("biopax-L2")) {
                String rdfFormat = getRdfFormatName(format);
                String mediaType = (format.contains("xml") || format.toLowerCase().startsWith("biopax"))? MediaType.APPLICATION_XML : MediaType.TEXT_PLAIN;

                psidev.psi.mi.xml.model.EntrySet entrySet = createEntrySet(query, firstResult, maxResults);

                StringWriter sw = new StringWriter();
                String output = "";

                PsimiRdfConverter rdfConverter = new PsimiRdfConverter();
                try {
                    rdfConverter.convert(entrySet, rdfFormat , sw);
                    output = sw.toString();
                } catch (Exception e) {
                    return formatNotSupportedResponse(format);
                }
                finally {
                    sw.close();
                    rdfConverter.close();
                }

                int count = 0;
                if (entrySet.getEntries() != null && !entrySet.getEntries().isEmpty()){
                    for (psidev.psi.mi.xml.model.Entry entry : entrySet.getEntries()){
                        if (entry.getInteractions() != null && !entry.getInteractions().isEmpty()){
                            count += entry.getInteractions().size();
                        }
                    }
                }
                return prepareResponse(Response.status(200).type(mediaType), output, count, isCompressed).build();

            } else {
                int count = count(query);
                
                if (RETURN_TYPE_COUNT.equalsIgnoreCase(format)) {
                    return prepareResponse(Response.status(200).type(MediaType.TEXT_PLAIN),
                            count, count, false)
                            .build();
                } else if (RETURN_TYPE_XGMML.equalsIgnoreCase(format)) {
                    XgmmlStreamingOutput xgmml = new XgmmlStreamingOutput(this.psicquicService, query, firstResult, maxResults);

                    String fixedQuery = query;
                    if (fixedQuery.contains("&")){
                        fixedQuery = query.substring(0, query.indexOf("&"));
                    }
                    fixedQuery = fixedQuery.replaceAll("q=", "");
                    fixedQuery = fixedQuery.replaceAll(":","_");
                    fixedQuery = fixedQuery.replaceAll(" ","_");
                    fixedQuery = fixedQuery.replaceAll("\\(","");
                    fixedQuery = fixedQuery.replaceAll("\\)","");

                    String name = fixedQuery.substring(0, Math.min(10, fixedQuery.length())) + ".xgmml";

                    return prepareResponse(Response.status(200).type(MediaType.APPLICATION_XML).header("Content-Disposition", "attachment; filename="+name),
                            xgmml, count, isCompressed)
                            .build();


                } else if (RETURN_TYPE_MITAB25.equalsIgnoreCase(format) || format == null) {
                    PsicquicStreamingOutput result = new PsicquicStreamingOutput(psicquicService, query, firstResult, maxResults, isCompressed);
                    return prepareResponse(Response.status(200).type(MediaType.TEXT_PLAIN), result,
                           result.countResults(), isCompressed).build();
                } else {
                    return formatNotSupportedResponse(format);
                }
            }
        } catch (Throwable e) {
            throw new PsicquicServiceException("Problem creating output", e);
        }


    }

    private Response formatNotSupportedResponse(String format) {
        return Response.status(406).type(MediaType.TEXT_PLAIN).entity(new GenericEntity<String>("Format not supported: " + format) {
        }).build();
    }

    private Response.ResponseBuilder prepareResponse(Response.ResponseBuilder responseBuilder, Object entity, long totalCount, boolean compressed) throws IOException {
        if (compressed) {
            if (entity instanceof InputStream) {
                CompressedStreamingOutput streamingOutput = new CompressedStreamingOutput((InputStream)entity);
                responseBuilder.entity(new GenericEntity<CompressedStreamingOutput>(streamingOutput){});
            } else if (entity instanceof String) {
                ByteArrayInputStream inputStream = new ByteArrayInputStream(((String)entity).getBytes());
                try{
                    CompressedStreamingOutput streamingOutput = new CompressedStreamingOutput(inputStream);
                    responseBuilder.entity(new GenericEntity<CompressedStreamingOutput>(streamingOutput){});
                }
                finally {
                    inputStream.close();
                }
            } else if (entity instanceof EntrySet) {
                ByteArrayOutputStream baos = new ByteArrayOutputStream();

                PsimiXmlWriter254 xmlWriter254 = new PsimiXmlWriter254();
                try {
                    xmlWriter254.marshall((EntrySet)entity, baos);

                    ByteArrayInputStream inputStream = new ByteArrayInputStream(baos.toByteArray());
                    try{
                        CompressedStreamingOutput streamingOutput = new CompressedStreamingOutput(inputStream);
                        responseBuilder.entity(new GenericEntity<CompressedStreamingOutput>(streamingOutput){});
                    }
                    finally {
                        inputStream.close();
                    }

                } catch (Throwable e) {
                    throw new IOException("Problem marshalling XML", e);
                }
                finally {
                    baos.close();
                }

            } else {
                responseBuilder.entity(new GenericEntity<Object>(entity){});
            }

            responseBuilder.header("Content-Encoding", "gzip");
        } else {
            responseBuilder.entity(new GenericEntity<Object>(entity){});
        }

        prepareHeaders(responseBuilder).header("X-PSICQUIC-Count", String.valueOf(totalCount));
        

        return responseBuilder;
    }
    
    public Response.ResponseBuilder prepareHeaders(Response.ResponseBuilder responseBuilder) {
        responseBuilder.header("X-PSICQUIC-Impl", config.getImplementationName());
        responseBuilder.header("X-PSICQUIC-Impl-Version", config.getVersion());
        responseBuilder.header("X-PSICQUIC-Spec-Version", config.getRestSpecVersion());
        responseBuilder.header("X-PSICQUIC-Supports-Compression", Boolean.TRUE);
        responseBuilder.header("X-PSICQUIC-Supports-Formats", StringUtils.join(SUPPORTED_REST_RETURN_TYPES, ", "));

        return responseBuilder;
    }

    private String getRdfFormatName(String format) {
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

    private psidev.psi.mi.xml.model.EntrySet createEntrySet(String query, int firstResult, int maxResults) throws ConverterException, PsicquicServiceException, NotSupportedMethodException, NotSupportedTypeException {
        EntrySetConverter converter = new EntrySetConverter();
        converter.setDAOFactory(new InMemoryDAOFactory());
        psidev.psi.mi.xml.model.EntrySet entrySet = converter.fromJaxb(getByQueryXml(query, firstResult, maxResults));
        return entrySet;
    }

    public Response getSupportedFormats() throws PsicquicServiceException, NotSupportedMethodException, NotSupportedTypeException {
        return Response.status(200)
                .type(MediaType.TEXT_PLAIN)
                .entity(new GenericEntity<String>(StringUtils.join(SUPPORTED_REST_RETURN_TYPES, "\n")){}).build();
    }

    public Object getProperty(String propertyName) {
        final String val = config.getProperties().get(propertyName);

        if (val == null) {
            return Response.status(404)
                .type(MediaType.TEXT_PLAIN)
                .entity(new GenericEntity<String>("Property not found: "+propertyName){}).build();
        }

         return Response.status(200)
                .type(MediaType.TEXT_PLAIN)
                .entity(new GenericEntity<String>(val){}).build();
    }

    public Response getProperties() {
        StringBuilder sb = new StringBuilder(256);

        for (Map.Entry entry : config.getProperties().entrySet()) {
            sb.append(entry.getKey()).append("=").append(entry.getValue()).append("\n");
        }

        return Response.status(200)
                .type(MediaType.TEXT_PLAIN)
                .entity(new GenericEntity<String>(sb.toString()){}).build();
    }

    public String getVersion() {
        return config.getVersion();
    }

    public psidev.psi.mi.xml254.jaxb.EntrySet getByQueryXml(String query,
                                  int firstResult,
                                  int maxResults) throws PsicquicServiceException, NotSupportedMethodException, NotSupportedTypeException {
        RequestInfo reqInfo = new RequestInfo();
        reqInfo.setResultType("psi-mi/xml25");

        try {
            reqInfo.setFirstResult(firstResult);
        } catch (NumberFormatException e) {
            throw new PsicquicServiceException("firstResult parameter is not a number: "+firstResult);
        }

        try {
            reqInfo.setBlockSize(maxResults);
        } catch (NumberFormatException e) {
            throw new PsicquicServiceException("maxResults parameter is not a number: "+maxResults);
        }

        QueryResponse response = psicquicService.getByQuery(query, reqInfo);

        return response.getResultSet().getEntrySet();
    }

    private int count(String query) throws NotSupportedTypeException, NotSupportedMethodException, PsicquicServiceException {
        RequestInfo reqInfo = new RequestInfo();
        reqInfo.setResultType("count");
        QueryResponse response = psicquicService.getByQuery(query, reqInfo);
        return response.getResultInfo().getTotalResults();
    }

    private String createQueryValue(String interactorAc, String db) {
        StringBuilder sb = new StringBuilder(256);
        if (db.length() > 0) sb.append('"').append(db).append(':');
        sb.append(interactorAc);
        if (db.length() > 0) sb.append('"');

        return sb.toString();
    }
}
