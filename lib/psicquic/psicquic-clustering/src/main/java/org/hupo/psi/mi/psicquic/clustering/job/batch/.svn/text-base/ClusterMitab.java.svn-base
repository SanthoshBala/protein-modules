package org.hupo.psi.mi.psicquic.clustering.job.batch;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.tab.model.builder.MitabWriterUtils;
import psidev.psi.mi.tab.model.builder.PsimiTabVersion;
import uk.ac.ebi.enfin.mi.cluster.ClusterContext;
import uk.ac.ebi.enfin.mi.cluster.Encore2Binary;
import uk.ac.ebi.enfin.mi.cluster.EncoreInteraction;
import uk.ac.ebi.enfin.mi.cluster.InteractionClusterAdv;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheStrategy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public class ClusterMitab {

    private static final Log log = LogFactory.getLog( ClusterMitab.class );

    public static void main( String[] args ) throws IOException {

        String miql = "";
        List<String> services = Arrays.asList( "IntAct", "Mint" );
        File dataLocationFile = new File("");

        final String evidenceMitabFile = args[0];
        final String clusteredMitabFile = args[1];

        ClusterContext.getInstance().setCacheStrategy( CacheStrategy.IN_MEMORY );

        InteractionClusterAdv iC = new InteractionClusterAdv();

        iC.addQueryAcc( miql );
        for ( String service : services ) {
            iC.addQuerySource( service );
        }

        // Run clustering
        iC.setMappingIdDbNames( "uniprotkb,irefindex,ddbj/embl/genbank,refseq,chebi" );
        iC.runService();

        Map<Integer, EncoreInteraction> interactionMapping = iC.getInteractionMapping();
        Encore2Binary iConverter = new Encore2Binary( iC.getMappingIdDbNames() );

        log.debug( "-------- CLUSTERED MITAB (" + interactionMapping.size() + ") --------" );


        String jobId = "";
        final File jobDirectory = new File( dataLocationFile, jobId );
        boolean created = jobDirectory.mkdirs();
        if ( !created ) {
            // TODO handle error
        }

        log.debug( "Using Job directory: " + jobDirectory.getAbsolutePath() );

        final String mitabFilename = jobId + ".tsv";
        final File mitabFile = new File( jobDirectory, mitabFilename );
        log.debug( "MITAB file: " + mitabFile.getAbsolutePath() );
        BufferedWriter out = new BufferedWriter( new FileWriter( mitabFile ) );

        for ( Map.Entry<Integer, EncoreInteraction> entry : interactionMapping.entrySet() ) {
            final EncoreInteraction ei = entry.getValue();
            final BinaryInteraction bi = iConverter.getBinaryInteraction( ei );
            final String mitab = MitabWriterUtils.buildLine(bi, PsimiTabVersion.v2_5);

            log.trace( mitab );
            out.write( mitab);
        }
    }
}
