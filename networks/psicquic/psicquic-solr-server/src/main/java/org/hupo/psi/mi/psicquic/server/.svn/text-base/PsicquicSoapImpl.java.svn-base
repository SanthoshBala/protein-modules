package org.hupo.psi.mi.psicquic.server;

/* =============================================================================
 # $Id::                                                                       $
 # Version: $Rev::                                                             $
 #==============================================================================
 #
 # PsicquicPortImpl: implementation of PSICQUIC 1.1 SOAP service 
 #
 #=========================================================================== */

import java.util.*;
import java.io.InputStream;

import javax.annotation.*;
import javax.jws.WebService;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.hupo.psi.mi.psq.*;
import org.hupo.psi.mi.psicquic.*;

import org.hupo.psi.mi.psicquic.server.index.*;
import org.hupo.psi.mi.psicquic.server.store.*;

import org.hupo.psi.mi.psicquic.util.JsonContext;

@WebService( name = "PsicquicService", 
             targetNamespace = "http://psi.hupo.org/mi/psicquic",
             serviceName = "PsicquicService",
             portName = "IndexBasedPsicquicServicePort",
             endpointInterface = "org.hupo.psi.mi.psq.PsqPort")
//             wsdlLocation = "/WEB-INF/wsdl/psicquic11.wsdl")

public class PsicquicSoapImpl implements PsqPort {
    
    org.hupo.psi.mi.psq.ObjectFactory psqOF =
        new org.hupo.psi.mi.psq.ObjectFactory();
    
    PsqContext psqContext;
    
    public void setPsqContext( PsqContext context ){
        psqContext = context;
    }


    PsicquicServer psqServer;

    public void setPsqServer( PsicquicServer server ){
        psqServer = server;
    }
    
    //--------------------------------------------------------------------------

    private void initialize() {
        initialize( false );
    }
    
    //--------------------------------------------------------------------------
    
    private void initialize( boolean force) {

	Log log = LogFactory.getLog( this.getClass() );
	log.info( " psqContext=" + psqContext );
    }

    //==========================================================================
    // WEB SERVICE OPERATIONS
    //=======================
   
    public QueryResponse getByQuery( String query,
                                     RequestInfo infoRequest )
        throws NotSupportedMethodException, 
               NotSupportedTypeException, 
               PsicquicServiceException {
        
        Log log = LogFactory.getLog( this.getClass() );
        log.info( "PsqPortImpl: getByQuery: context =" + psqContext);
        log.info( "PsqPortImpl: getByQuery: q=" + query );
        
        long firstResult = -1;
        long blockSize = -1;
        String viewType = null; 
        
        if( infoRequest != null ){              
            log.info( "PsqPortImpl: FR=" + infoRequest.getFirstResult() 
                      + " BS=" + infoRequest.getBlockSize() );
            
            firstResult = infoRequest.getFirstResult();           
            blockSize = infoRequest.getBlockSize();
               
            if( infoRequest.getResultType() != null ){
                viewType = infoRequest.getResultType();   
            } 
        }
        
        ResultSet qrs = psqServer.getByQuery( query, viewType,
                                              firstResult, blockSize );
        
        QueryResponse qr = psqOF.createQueryResponse();
        qr.setResultSet( psqOF.createResultSet() );
        
        String mitab =  psqServer.getHeader( viewType );
        
	log.info( "getByQuery: rs="+ qrs); 
        
        for( Iterator i = qrs.getResultList().iterator(); i.hasNext(); ){
            String record = (String) i.next();
            mitab += record + "\n";
        }
        mitab += psqServer.getFooter( viewType );

        qr.getResultSet().setMitab( mitab );     

        return qr;
        
    };

    //--------------------------------------------------------------------------
    
    public QueryResponse getByInteractor( DbRef dbRef,
                                          RequestInfo infoRequest )
        throws NotSupportedMethodException, 
               NotSupportedTypeException, 
               PsicquicServiceException {

        String query = psqServer.buildQuery( "identifier", 
                                             dbRef.getDbAc(), dbRef.getId() );        
        return getByQuery( query, infoRequest );
    };
    
    //--------------------------------------------------------------------------

    public QueryResponse getByInteractorList( List<DbRef> dbRef,
                                              RequestInfo infoRequest,
                                              String operand )
        throws NotSupportedMethodException, 
               NotSupportedTypeException, 
               PsicquicServiceException {

        List<String> dbl = new ArrayList();
        List<String> acl = new ArrayList();

        for( Iterator<DbRef> i = dbRef.iterator(); i.hasNext(); ){
            DbRef cref = i.next();
            dbl.add( cref.getDbAc());
            acl.add( cref.getId());
        }

        String query = psqServer.buildQuery( "identifier", dbl, acl, operand );
        return getByQuery( query, infoRequest );
    };

    //--------------------------------------------------------------------------
    
    public QueryResponse getByInteraction( DbRef dbRef,
                                           RequestInfo infoRequest )
        throws NotSupportedMethodException, 
               NotSupportedTypeException, 
               PsicquicServiceException {
        
        String query = psqServer.buildQuery( "interaction_id",
                                             dbRef.getDbAc(), dbRef.getId() );
        return getByQuery( query, infoRequest );        
    };

    //--------------------------------------------------------------------------
    
    public QueryResponse getByInteractionList( List<DbRef> dbRef,
                                               RequestInfo infoRequest )
        throws NotSupportedMethodException, 
               NotSupportedTypeException, 
               PsicquicServiceException{
       
        List<String> dbl = new ArrayList();
        List<String> acl = new ArrayList();

        for( Iterator<DbRef> i = dbRef.iterator(); i.hasNext(); ){
            DbRef cref = i.next();
            dbl.add( cref.getDbAc());
            acl.add( cref.getId());
        }
        
        String query = psqServer.buildQuery( "identifier", dbl, acl, "OR" );
        return getByQuery( query, infoRequest );
    };

    //==========================================================================
    // META DATA
    //==========

    public List<String> getSupportedReturnTypes(){
        return psqServer.getSupportedReturnTypes( "soap" );
    };
    
    //--------------------------------------------------------------------------
    
    public String getVersion(){
        return psqServer.getVersion( "soap" );
    }

    //--------------------------------------------------------------------------
    
    public List<String> getSupportedDbAcs(){
        return psqServer.getSupportedDbAcs( "soap" );
    }
    
    //--------------------------------------------------------------------------

    public String getProperty( String property ){     
        return psqServer.getProperty( "soap", property );
    }
    
    //--------------------------------------------------------------------------
    
    public List<Property> getProperties(){
        
        Set<Map.Entry> props = psqServer.getProperties( "soap" );
        List<Property> pl = new ArrayList<Property>();
        
        for( Iterator<Map.Entry> pi = props.iterator(); pi.hasNext(); ){
            Map.Entry me = pi.next();
            
            Property p = psqOF.createProperty();
        
            p.setKey( (String) me.getKey() );
            p.setValue( (String) me.getValue() );
            pl.add( p );
        }
        return pl;
    }
}
