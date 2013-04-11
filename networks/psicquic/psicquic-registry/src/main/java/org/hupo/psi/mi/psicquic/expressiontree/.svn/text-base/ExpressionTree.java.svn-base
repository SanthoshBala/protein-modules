package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

import java.util.List;
import java.util.Set;

/**
 * 
 */

/**
 * @author Erik Pfeiffenberger
 *
 */

public class ExpressionTree {
	
	private ExpressionNode rootNode;
	private List<String> tags;
	private String expression;
	

	private SelfDiscoveringOntologyTree miOntologyTree;
	private boolean lazyBuild;
	
	private ExpressionScanner scanner;
	
	
	public ExpressionTree(String expression, SelfDiscoveringOntologyTree miOntologyTree, boolean lazyBuild) throws ParseExpressionException{
		this.expression = expression;
		this.lazyBuild = lazyBuild;
		this.miOntologyTree = miOntologyTree;
		
		if(!this.lazyBuild){
			rootNode = parse();
		} else {
			rootNode = null;
		}
	}
	
	
	private ExpressionNode parse() throws ParseExpressionException {
		ExpressionNode rootNode = null;
		if(isCommaSeparatedList()){
			rootNode = parseCommaseperatedList();
		} else {
			scanner = new ExpressionScanner(this.expression);
			rootNode = parseExpression();
		}
		
		return rootNode;
	}

	private ExpressionNode parseCommaseperatedList() throws ParseExpressionException {
		String [] elements = expression.split(",");
		String newExpression = "(" + elements[0].trim();
		
		for (int i=1; i< elements.length; i++){
			newExpression = newExpression + " or " + elements[i].trim();
		}
		
		newExpression = newExpression +")";
		
		expression = newExpression;
		scanner = new ExpressionScanner(this.expression);
		return parseExpression();
	}

	public boolean isCommaSeparatedList() {
		//comma separated list, e.g. "MI:0000 , MI:0001, some free text"
		return expression.matches("(.+,.+)+");
	
	}

	private ExpressionNode parseExpression() throws ParseExpressionException {
		ExpressionNode rootNode = null;

		rootNode = parseTerm();
		scanner.nextSymbol();	
		while(isLogicalOp()){
			if(scanner.isAnd()){
				ExpressionNode oldRootNode = rootNode;
				rootNode = new ExpressionAndNode(oldRootNode, parseTerm());
			} else if(scanner.isOr()){
				ExpressionNode oldRootNode = rootNode;
				rootNode = new ExpressionOrNode(oldRootNode, parseTerm());
			} else {
					throw new ParseExpressionException();
			}
		}
		
		
		return rootNode;
	}

	private boolean isLogicalOp() {
		return (scanner.isAnd() || scanner.isOr());
	}

	private ExpressionNode parseTerm() throws ParseExpressionException {
		ExpressionNode rootNode = null;
		
		scanner.nextSymbol();
		while(scanner.isNot()){
			
			ExpressionNode oldRoot = rootNode;
			rootNode = new ExpressionNotNode(oldRoot);
			scanner.nextSymbol();
		}
		
		if(scanner.isMiIdentifier()){
			if(rootNode != null){
				ExpressionNotNode tmpNode = (ExpressionNotNode)rootNode;
				while(tmpNode.getChild() != null){
					tmpNode = (ExpressionNotNode)tmpNode.getChild();
				}
				tmpNode.setChild(new ExpressionMiNode(scanner.getCurrSymbol()));
			} else {
				rootNode = new ExpressionMiNode(scanner.getCurrSymbol());
			}
		} else if (scanner.isFreeText()){
			if(rootNode != null){
				ExpressionNotNode tmpNode = (ExpressionNotNode)rootNode;
				while(tmpNode.getChild() != null){
					tmpNode = (ExpressionNotNode)tmpNode.getChild();
				}
				tmpNode.setChild(parseFreeText());
				
			} else {
				rootNode = parseFreeText();
			}
		} else if (scanner.isLeftPar()){
			if(rootNode != null){
				ExpressionNotNode tmpNode = (ExpressionNotNode)rootNode;
				while(tmpNode.getChild() != null){
					tmpNode = (ExpressionNotNode)tmpNode.getChild();
				}
				tmpNode.setChild(parseExpression());
			} else {
				rootNode = parseExpression();
			}
			
			scanner.nextSymbol();
			if(!scanner.isRightPar()){
				throw new ParseExpressionException();
			}
		} else {
			throw new ParseExpressionException();
		}
		
		
		return rootNode;
	}

	private ExpressionNode parseFreeText() throws ParseExpressionException {
		Set<String>ids = miOntologyTree.getIdByName(scanner.getCurrSymbol());
		
		ExpressionNode rootNode = null;
		
		if(ids.size() > 0){
			Object [] idsArr = ids.toArray();
			String tmpExpression = (String)idsArr[0];
			
			for (int i = 1; i< idsArr.length; ++i){
				tmpExpression = tmpExpression + " or " + (String)idsArr[i];
			}
			
			ExpressionTree tmpTree = new ExpressionTree(tmpExpression, miOntologyTree, false);
			
			rootNode = tmpTree.getRootNode();
			
		} else {
            // TODO FIXME if I give a tag that doesn't have an MI identifier, we return all services.
            // TODO collect these non matching terms and print a warning message on the resulting HTML page.
			rootNode = new ExpressionTrueNode();
		}
		
		return rootNode;
	}


	public boolean evaluate (ServiceType service) throws ParseExpressionException {
		if (rootNode == null){
			rootNode = parse();
		}
		
		return rootNode.evaluate(service, miOntologyTree);
	}
	
	public ExpressionNode getRootNode() throws ParseExpressionException{		
		if (rootNode == null){
			rootNode = parse();
		}
		
		return rootNode;
	}
	
	public void printTree(){
		doPrintTree(rootNode, 0);
	}


	private void doPrintTree(ExpressionNode node, int depth) {
		if(node != null){
			if(node.getLeftChild() != null){
				doPrintTree(node.getLeftChild(), depth+1);
			}
			
			if(node.getRightChild() != null){
				doPrintTree(node.getRightChild(), depth+1);
			}
		}
		
	}
	
	private String repeatStr(String str, int nTimes){
		String returnVal = "";
		
		if(nTimes >= 0){
			for(int i=0; i<=nTimes; i++){
				returnVal = returnVal+str;
			}
		} else {
			returnVal = str;
		}
		
		return returnVal;
	}

}
