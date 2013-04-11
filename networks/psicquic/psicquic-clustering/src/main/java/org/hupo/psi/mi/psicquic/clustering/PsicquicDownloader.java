package org.hupo.psi.mi.psicquic.clustering;

import org.hupo.psi.mi.psicquic.registry.ServiceType;
import org.hupo.psi.mi.psicquic.registry.client.PsicquicRegistryClientException;
import org.hupo.psi.mi.psicquic.registry.client.registry.DefaultPsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.registry.client.registry.PsicquicRegistryClient;
import org.hupo.psi.mi.psicquic.wsclient.PsicquicSimpleClient;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;


/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class PsicquicDownloader {


    public static Collection<Service> getServiceURLs( List<String> serviceNames ) throws PsicquicRegistryClientException {
        PsicquicRegistryClient registry = new DefaultPsicquicRegistryClient();
        final List<ServiceType> services = registry.listActiveServices();

        System.out.println( "Found " + services.size() + " active service(s)." );
        Collection selectedServices = new ArrayList<Service>( services.size() );
        for ( ServiceType st : services ) {
            final String name = st.getName();
            if ( serviceNames.contains( name ) ) {
                Service s = new Service( name );
                s.setRestUrl( st.getRestUrl() );
                selectedServices.add( s );

                System.out.println( "Added Service '" + name + "' with REST URL: '" + st.getRestUrl() + "'" );
            }
        }

        return selectedServices;
    }

    public static void main( String[] args ) throws Exception {
        // get a REST URl from the registry http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS

        final Collection<Service> services = getServiceURLs( Arrays.asList( "IntAct", "MINT", "foobar" ) );

        final String query = "brca2";

        for ( Service service : services ) {

            PsicquicSimpleClient client = new PsicquicSimpleClient( service.getRestUrl() );

            // count the results
            final long count = client.countByQuery( query );

            System.out.println( "\n\n" + service.getName() + " (" + count + " interactions)" );
            System.out.println( "-----------------------------------------" );

            try {
                final InputStream result = client.getByQuery( query, PsicquicSimpleClient.MITAB25 );

                BufferedReader in = new BufferedReader( new InputStreamReader( result ) );

                String line;

                while ( ( line = in.readLine() ) != null ) {
                    System.out.println( line );
                }

                in.close();
                result.close();
            } catch ( IOException e ) {
                e.printStackTrace();
            }

        } // services

    } // main
}
