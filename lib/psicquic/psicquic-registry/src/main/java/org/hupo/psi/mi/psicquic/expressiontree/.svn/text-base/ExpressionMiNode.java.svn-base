package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

import java.util.List;
/**
 * 
 */

/**
 * @author Erik Pfeiffenberger
 *
 */
public class ExpressionMiNode implements ExpressionNode{
	
	private String miId;
	
	public ExpressionMiNode(String miId) {
	
		this.miId = miId.toLowerCase();
	}

	public boolean evaluate(ServiceType service, SelfDiscoveringOntologyTree sdoTree) {
		boolean result = false;
		List<String> serviceTags = service.getTags();
		
		for(String sTag : serviceTags){
			if(sTag.equalsIgnoreCase(miId)){
				result = true;
				break;
			}
		}
		
		if(!result){
			for(String sTag : serviceTags){
				if(sdoTree.isChildTag(miId, sTag)){
					result = true;
					break;
				}
			}
		}
		
		return result;
	}

	public String getRepresentation() {
		return miId;
	}

	public ExpressionNode getLeftChild() {
		return null;
	}

	public ExpressionNode getRightChild() {
		return null;
	}



}
