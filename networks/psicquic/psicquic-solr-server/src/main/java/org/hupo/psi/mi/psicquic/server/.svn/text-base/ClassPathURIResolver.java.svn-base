package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # Classpath aware URIResolver:
 #
 #     Attempts to resolve classpath: scheme to current classpath, if fails
 #     falls back to filesystem (after dropping scheme prefix). 
 #     Scheme-less URIs are firsr resolved within the filesystem and then, 
 #     in case of failure, against current classpath
 #
 #=========================================================================== */

//import java.util.*;
//import java.net.*;
import java.io.*;

import javax.xml.transform.Source;
import javax.xml.transform.stream.StreamSource;
import javax.xml.transform.URIResolver;
import javax.xml.transform.TransformerException;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class ClassPathURIResolver implements URIResolver{
    
    private URIResolver defURIResolver = null;
    
    public ClassPathURIResolver( URIResolver resolver ){
        defURIResolver = resolver;
    }
    
    public Source resolve( String href,  String base)
        throws TransformerException {
        
        if( href == null ){ 
            throw new TransformerException( "URIResolver: href missing" );
        }
        
        if( href.startsWith( "classpath:" ) ){
            href = href.substring( 10 ); 
            
            InputStream is =  this.getClass().getClassLoader()
                .getResourceAsStream( href );
            
            if( is != null ){
                return new StreamSource( is );
            }   
        } else {
            try{
                InputStream fis = new FileInputStream( href );
                return new StreamSource( fis );
            } catch(FileNotFoundException fx ){
                InputStream is =  this.getClass().getClassLoader()
                    .getResourceAsStream( href );

                if( is != null ){
                    return new StreamSource( is );
                }
            }
        }
        return defURIResolver.resolve( href, base );
    }
}
