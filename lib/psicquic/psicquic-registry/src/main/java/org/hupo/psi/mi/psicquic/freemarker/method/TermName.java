/**
 * 
 */
package org.hupo.psi.mi.psicquic.freemarker.method;

import freemarker.template.TemplateMethodModel;
import freemarker.template.TemplateModelException;
import org.hupo.psi.mi.psicquic.ols.client.OntologyTerm;
import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;

import java.util.List;

/**
 * @author Erik Pfeiffenberger
 *
 */

public class TermName implements TemplateMethodModel {

	SelfDiscoveringOntologyTree miOntologyTree;
	
	public SelfDiscoveringOntologyTree getMiOntologyTree() {
		return miOntologyTree;
	}

	public void setMiOntologyTree(SelfDiscoveringOntologyTree miOntologyTree) {
		this.miOntologyTree = miOntologyTree;
	}

	public TermName(SelfDiscoveringOntologyTree miOntologyTree) {
		this.miOntologyTree = miOntologyTree;
	}
	
	public Object exec(List arg0) throws TemplateModelException {
		String termName;

		if(arg0.size() != 1){
			throw new TemplateModelException("Wrong arguments");

		} else {
			String id = (String)arg0.get(0);
			
			OntologyTerm term = miOntologyTree.getTerm(id);
			
			termName = term.getName();
		}
		return termName;
	}

}
