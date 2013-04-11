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

import org.hupo.psi.mi.psicquic.DbRef;
import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.RequestInfo;
import org.hupo.psi.mi.psicquic.ws.config.PsicquicConfig;
import org.junit.*;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import psidev.psi.mi.search.Searcher;

import java.io.File;
import java.io.InputStream;
import java.util.Arrays;

/**
 * TODO comment that class header
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id$
 */
public class IndexBasedPsicquicServiceTest {

    private static IndexBasedPsicquicService service;

    @BeforeClass
    public static void beforeClass() throws Exception {
        InputStream mitabStream = IndexBasedPsicquicServiceTest.class.getResourceAsStream("/META-INF/brca2.mitab.txt");
        File indexDir = new File("target", "brca-mitab.index");

        Searcher.buildIndex(indexDir, mitabStream, true, true);

        ClassPathXmlApplicationContext context = new ClassPathXmlApplicationContext(new String[] {"/META-INF/beans.spring.test.xml"});
        PsicquicConfig config = (PsicquicConfig)context.getBean("psicquicConfig");
        config.setIndexDirectory(indexDir.toString());

	    service = (IndexBasedPsicquicService)context.getBean("indexBasedPsicquicService");
    }

    @AfterClass
    public static void afterClass() throws Exception {
        service = null;
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testGetByInteractor() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        DbRef dbRef = new DbRef();
        dbRef.setId("FANCD1");

        final QueryResponse response = service.getByInteractor(dbRef, info);

        Assert.assertEquals(12, response.getResultInfo().getTotalResults());
        Assert.assertEquals(12, response.getResultSet().getMitab().split("\n").length);
    }

    @Test
    public void testGetByInteraction() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        DbRef dbRef = new DbRef();
        dbRef.setId("EBI-1372935");

        final QueryResponse response = service.getByInteraction(dbRef, info);

        Assert.assertEquals(2, response.getResultInfo().getTotalResults());
    }

    @Test
    public void testGetByInteractorList_operandOR() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        DbRef dbRef1 = new DbRef();
        dbRef1.setId("BRCA2");
        DbRef dbRef2 = new DbRef();
        dbRef2.setId("FANCD2");

        final QueryResponse response = service.getByInteractorList(Arrays.asList(dbRef1, dbRef2), info, "OR");

        Assert.assertEquals(14, response.getResultInfo().getTotalResults());

    }

    @Test
    public void testGetByInteractorList_operandAND() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        DbRef dbRef1 = new DbRef();
        dbRef1.setId("BRCA2");
        DbRef dbRef2 = new DbRef();
        dbRef2.setId("FANCD2");

        final QueryResponse response = service.getByInteractorList(Arrays.asList(dbRef1, dbRef2), info, "AND");

        Assert.assertEquals(2, response.getResultInfo().getTotalResults());

    }

    @Test
    public void testGetByInteractionList() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        DbRef dbRef1 = new DbRef();
        dbRef1.setId("EBI-1372935");
        DbRef dbRef2 = new DbRef();
        dbRef2.setId("EBI-1372896");

        final QueryResponse response = service.getByInteractionList(Arrays.asList(dbRef1, dbRef2), info);

        Assert.assertEquals(2, response.getResultInfo().getTotalResults());
    }

    @Test
    public void testGetByQuery() throws Exception {
        RequestInfo info = new RequestInfo();
        info.setResultType("psi-mi/tab25");
        info.setBlockSize(50);

        final QueryResponse response = service.getByQuery("FANCD1", info);

        Assert.assertEquals(12, response.getResultInfo().getTotalResults());
        Assert.assertEquals(12, response.getResultSet().getMitab().split("\n").length);
    }

    @Test
    public void testGetVersion() {
        Assert.assertEquals("TEST.VERSION", service.getVersion());
    }

}
