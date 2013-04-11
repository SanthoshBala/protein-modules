package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # XsltTransformer: transform the incoming xml file into solr documents
 #
 #=========================================================================== */

import java.util.*;
import java.net.*;
import java.io.*;

import org.w3c.dom.*;
import javax.xml.parsers.*;

import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.StreamSource;
import javax.xml.bind.util.JAXBResult;

import org.apache.solr.common.SolrInputDocument;

public class XsltTransformer implements PsqTransformer{
    
    private Transformer psqtr = null;

    private String fileName = "";
    private NodeList domList = null;
    private int domPos = -1;
    
    public XsltTransformer( String xslt ){
	
	try {
	    DocumentBuilderFactory
		dbf = DocumentBuilderFactory.newInstance();
	    dbf.setNamespaceAware( true );
	    
	    DocumentBuilder db = dbf.newDocumentBuilder();
	    
	    File xslFile = new File( xslt );
	    InputStream xslIStr = null;
	    
	    if( !xslFile.canRead() ){
		xslIStr = this.getClass().getClassLoader()
		    .getResourceAsStream( xslt );
	    } else {
		xslIStr = new FileInputStream( xslt );
	    }
	    
	    Document xslDoc = db.parse( xslIStr );
	    
	    DOMSource xslDomSource = new DOMSource( xslDoc );
	    TransformerFactory
		tFactory = TransformerFactory.newInstance();
	    
	    psqtr = tFactory.newTransformer( xslDomSource );
            
	} catch( Exception ex ) {
	    ex.printStackTrace();
	}
    }

    public void start( String fileName ){
	this.fileName = fileName;	
    }
    
    public void start( String fileName, InputStream is ){
	
	try {
	    this.fileName = fileName;
	    
            StreamSource ssNative = new StreamSource( is );
            DOMResult domResult = new DOMResult();
	    
            //transform into dom
            //------------------
            
            psqtr.clearParameters();
            psqtr.setParameter( "file", fileName );
            psqtr.transform( ssNative, domResult );
            
	    Node domNode = domResult.getNode();
	    
	    if( domNode != null ){
	        domList = domNode.getFirstChild().getChildNodes();
	    }    
	    
	} catch( Exception ex ) {
	    ex.printStackTrace();
	}
    }

    public boolean hasNext(){
	if( domList == null | domPos >= domList.getLength() ) return false;
	return true;
    }
    
    public SolrInputDocument next(){
	
	if( domList != null ){
	    domPos++;
	    
	    if( domPos < domList.getLength() ){
		
                if( domList.item( domPos ).getNodeName().equals("doc") ){
		    
                    String pid = fileName;

                    SolrInputDocument doc =  new SolrInputDocument();
                    NodeList field = domList.item( domPos ).getChildNodes();
                    for( int j = 0; j< field.getLength() ;j++ ){
                        if( field.item(j).getNodeName().equals("field") ){
                            Element fe = (Element) field.item(j);
                            String name = fe.getAttribute("name");
                            String value = fe.getFirstChild().getNodeValue();
                            doc.addField( name,value );
			    
                            if( name.equals( "pid" ) ){
                                pid = value;
                            }
                        }
                    }	    			      
		    return doc;
		}
	    }                                    
	}
       	return null;
    }
}