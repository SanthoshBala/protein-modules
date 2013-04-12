package org.hupo.psi.mi.psicquic.stats;
import org.hupo.psi.mi.psicquic.stats.PsicquicService;
import org.hupo.psi.mi.psicquic.stats.PsicquicStatsCollector;
import org.junit.Assert;
import org.junit.Test;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Test data collection. Count interaction for all the psicquic services
 * and retrieve publication information for one sservice. We chose "MatrixDB"
 * since it is one of the smallest sources and it is present in both PSICQUIC
 * and IMEx.
 *
 * @author Rafael Jimenez (rafael@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class TestPsicquicStatsCollector {
    @Test
    public void testData(){
        try {
            List<String> sourcesToQuery = new ArrayList<String>();
            String sourceTest = "MatrixDB";
            sourcesToQuery.add(sourceTest);

            PsicquicStatsCollector collector = new PsicquicStatsCollector();

            List<PsicquicService> psicquicServices = collector.collectPsicquicServiceNames();

            Map<String, Long> db2interactionCount = collector.updatePsicquicInteractionsStats(psicquicServices); //data
            Map<String, Long> db2publicationsCount = collector.collectPsicquicPublicationsStats(psicquicServices, sourcesToQuery); //data

//            System.out.println("db2interactionCount MatrixDB: " + db2interactionCount.get(sourceTest));

            Assert.assertTrue(db2interactionCount.containsKey(sourceTest));
            Assert.assertTrue(db2interactionCount.get(sourceTest) > 220 );
            Assert.assertTrue(db2publicationsCount.containsKey(sourceTest));
            Assert.assertTrue(db2publicationsCount.get(sourceTest) > 35 );
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
