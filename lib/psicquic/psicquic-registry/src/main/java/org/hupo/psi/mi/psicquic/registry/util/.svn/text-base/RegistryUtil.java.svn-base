/**
 * 
 */
package org.hupo.psi.mi.psicquic.registry.util;

import org.hupo.psi.mi.psicquic.registry.Registry;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.InputStream;


/**
 * @author Erik Pfeiffenberger
 *
 */
public class RegistryUtil {

    public static Registry readRegistry(InputStream is) {
        Registry objRegistry;
        try {
            objRegistry = (Registry) readRegistryXml(is);
        } catch (Throwable e) {
            throw new IllegalStateException("Problem reading registry", e);
        }

        return objRegistry;
    }

    private static Object readRegistryXml(InputStream is) throws JAXBException {
        JAXBContext jc = JAXBContext.newInstance(Registry.class.getPackage().getName());
        Unmarshaller unmarshaller = jc.createUnmarshaller();
        return unmarshaller.unmarshal(is);
    }
}
