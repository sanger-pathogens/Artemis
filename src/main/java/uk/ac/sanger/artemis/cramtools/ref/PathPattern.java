/*******************************************************************************
 * Copyright 2013 EMBL-EBI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package uk.ac.sanger.artemis.cramtools.ref;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Formats a string pattern with an id. The following notations are support:
 * <p>
 * "%[precision]s"
 * <p>
 * "%s", expected at the end of the pattern or the last formatting notation.
 * <p>
 * The id supplied in the argument will be trimmed from the left so that the
 * last %s notation will receive the remainder of the id.
 * <p>
 * For example: A pattern "anything/%2s/%3s/%s" formatted for id 0123456789 will
 * result in "anything/01/234/56789"
 * 
 * @author vadim
 * 
 */
class PathPattern {

	private String pathPatternFormat;
	private Pattern pattern = Pattern.compile("%(\\d+)s");

	public PathPattern(String format) {
		this.pathPatternFormat = format;
	}

	public String format(String id) {
		Matcher matcher = pattern.matcher(pathPatternFormat);
		StringBuffer sb = new StringBuffer();
		int prevEnd = 0;
		while (matcher.find()) {
			String group = matcher.group();
			int number = Integer.valueOf(group.substring(1, group.length() - 1));
			sb.append(pathPatternFormat.substring(prevEnd, matcher.start()));
			sb.append(id.substring(0, number));
			id = id.substring(number);
			prevEnd = matcher.end();
		}
		sb.append(pathPatternFormat.substring(prevEnd, pathPatternFormat.length()));
		if (!sb.toString().contains("%s"))
			throw new RuntimeException("Invalid reference location pattern: " + pathPatternFormat);
		return sb.toString().replaceAll("%s", id);
	}

}
