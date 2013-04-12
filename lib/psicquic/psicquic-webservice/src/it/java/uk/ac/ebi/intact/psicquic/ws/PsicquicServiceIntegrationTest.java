package uk.ac.ebi.intact.psicquic.ws;

import org.apache.cxf.frontend.ClientProxyFactoryBean;
import org.apache.cxf.jaxws.JaxWsProxyFactoryBean;
import org.hupo.psi.mi.psicquic.DbRef;
import org.hupo.psi.mi.psicquic.PsicquicService;
import org.hupo.psi.mi.psicquic.RequestInfo;
import org.junit.Ignore;
import org.junit.Test;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import psidev.psi.mi.xml254.jaxb.EntrySet;


public class PsicquicServiceIntegrationTest {

	@Test @Ignore
	public void getVersion() throws Exception {
		ClientProxyFactoryBean factory = new JaxWsProxyFactoryBean();
		factory.setServiceClass(PsicquicService.class);
		factory.setAddress("http://localhost:8080/psicquic-ws/webservices/psicquic");
		PsicquicService client = (PsicquicService) factory.create();

		String response = client.getVersion();
		System.out.println("Received response from webservice: " + response);
	}

    @Test @Ignore
	public void getVersion2() throws Exception {
		//lookup client
	    ClassPathXmlApplicationContext context = new ClassPathXmlApplicationContext(new String[] {"spring-test.xml"});

	    PsicquicService client = (PsicquicService)context.getBean("psicquicServiceClient");

        System.out.println("Version: " +client.getVersion());

        DbRef dbRef = new DbRef();
        dbRef.setId("P1234");

        RequestInfo request = new RequestInfo();
        request.setBlockSize(50);
        request.setFirstResult(0);

        final EntrySet entrySet = client.getByInteractor(dbRef, request).getResultSet().getEntrySet();
        System.out.println("Version: "+entrySet.getLevel()+"."+entrySet.getVersion()+"."+entrySet.getMinorVersion());
        System.out.println("Entries: "+entrySet.getEntries().size());
    }

}