package org.gmod.schema.utils;

/**
 * Class to store terms, eg products, and how many occurences there are 
 * 
 * @author cp2
 */
public class CountedName {
	
	/**
	 * 
	 */
	private String name;
	
	/**
	 * 
	 */
	private int count;

	
	
	public CountedName(String name, int count) {
		this.name = name;
		this.count = count;
	}

	/**
	 * @return
	 */
	public int getCount() {
		return count;
	}

	/**
	 * @param count
	 */
	public void setCount(int count) {
		this.count = count;
	}

	/**
	 * @return
	 */
	public String getName() {
		return name;
	}

	/**
	 * @param name
	 */
	public void setName(String name) {
		this.name = name;
	}

}
