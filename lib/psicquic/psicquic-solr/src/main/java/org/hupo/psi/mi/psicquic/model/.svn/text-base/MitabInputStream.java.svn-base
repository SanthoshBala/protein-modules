package org.hupo.psi.mi.psicquic.model;

import org.apache.commons.io.IOUtils;
import org.apache.solr.common.SolrDocument;
import org.apache.solr.common.SolrDocumentList;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.Iterator;

/**
 * Input stream for MITAB
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>30/07/12</pre>
 */

public class MitabInputStream extends InputStream{

    private Iterator<SolrDocument> solrResultsIterator;
    private StringBuffer mitabLineBuffer;
    private String [] fieldNames;
    private InputStream mitabLineInputStream;

    private static final String NEW_LINE = System.getProperty("line.separator");
    private static final String COLUMN_SEPARATOR = "\t";
    private static final String FIELD_SEPARATOR = "|";
    private static final String FIELD_EMPTY = "-";

    public MitabInputStream(SolrDocumentList solrResults, String [] fieldNames){

        if (solrResults != null){
            solrResultsIterator = solrResults.iterator();
        }

        this.fieldNames = fieldNames;
        this.mitabLineBuffer = new StringBuffer();
    }
    @Override
    public int read() throws IOException {

        // initializes mitabLineInputStream if not done yet
        if (mitabLineInputStream == null){
            mitabLineInputStream = readNextLine();
        }

        // we reached the end of results
        if (mitabLineInputStream == null){
            return -1;
        }

        int charRead = mitabLineInputStream.read();

        // we need to read the next MITAB line because we reached the end of the current one
        while (charRead == -1 && mitabLineInputStream != null){
            mitabLineInputStream = readNextLine();

            if (mitabLineInputStream != null){
                charRead = mitabLineInputStream.read();
            }
        }

        return charRead;
    }

    @Override
    public int read(byte[] bytes) throws java.io.IOException {

        return read(bytes, 0, bytes.length);
    }

    @Override
    public int read(byte[] bytes, int off, int len) throws java.io.IOException {
        if (bytes == null){
            throw new NullPointerException("The array of bytes to read for the MitabInputStream is null.");
        }
        if (off < 0){
            throw new IndexOutOfBoundsException("The off value to read the mitab inputStream cannot be negative " + off);
        }
        if (len < 0){
            throw new IndexOutOfBoundsException("The len value to read the mitab inputStream cannot be negative " + len);
        }
        if (len < bytes.length - off ){
            throw new IndexOutOfBoundsException("The len value to read the mitab inputStream cannot be inferior to bytes.length - off. Length = " + len + ", bytes.length - off = " + Integer.toString(bytes.length - off));
        }
        if (len == 0){
            return 0;
        }

        // initializes mitabLineInputStream if not done yet
        if (mitabLineInputStream == null){
            mitabLineInputStream = readNextLine();
        }

        // we reached the end of results
        if (mitabLineInputStream == null){
            return -1;
        }

        int numberCharRead = Math.max(0, mitabLineInputStream.read(bytes, off, len));

        // we need to read the next MITAB line because we reached the end of the current one
        while (numberCharRead < len && mitabLineInputStream != null){
            mitabLineInputStream = readNextLine();

            if (mitabLineInputStream != null){
                int newNumberCharRead = Math.max(0,mitabLineInputStream.read(bytes, off + numberCharRead, len - numberCharRead));
                numberCharRead+=newNumberCharRead;
            }
        }

        if (numberCharRead == 0){
           return -1;
        }

        return numberCharRead;
    }

    @Override
    public int available() throws java.io.IOException {

        if (mitabLineInputStream != null) {
            return mitabLineInputStream.available();
        }

        return 0;
    }

    private InputStream readNextLine() throws IOException {

        if (solrResultsIterator != null && solrResultsIterator.hasNext() && fieldNames != null && fieldNames.length > 0){
            SolrDocument solrDoc = solrResultsIterator.next();

            // clear the mitabLineBuffer
            mitabLineBuffer.setLength(0);

            // close previous inputStream
            if (mitabLineInputStream != null){
                mitabLineInputStream.close();
            }

            int size = fieldNames.length;
            int index = 0;
            // one field name is one column
            for (String fieldName : fieldNames){
                // only one value is expected because it should not be multivalued
                Collection<Object> fieldValues = solrDoc.getFieldValues(fieldName);

                if (fieldValues == null || fieldValues.isEmpty()){
                    mitabLineBuffer.append(FIELD_EMPTY);
                }
                else {
                    Iterator<Object> valueIterator = fieldValues.iterator();
                    while (valueIterator.hasNext()){
                        mitabLineBuffer.append(String.valueOf(valueIterator.next()));

                        if (valueIterator.hasNext()){
                            mitabLineBuffer.append(FIELD_SEPARATOR);
                        }
                    }
                }

                index++;

                if (index < size){
                    mitabLineBuffer.append(COLUMN_SEPARATOR);
                }
            }

            mitabLineBuffer.append(NEW_LINE);

            return IOUtils.toInputStream(mitabLineBuffer);
        }

        return null;
    }

    @Override
    public void close() throws IOException {
        mitabLineBuffer.setLength(0);

        if (mitabLineInputStream != null){
            mitabLineInputStream.close();
        }
        mitabLineInputStream = null;
    }
}
