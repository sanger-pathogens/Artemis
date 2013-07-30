/* DatabaseLocationParser.java
 *
 * created: Jun 2013
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */
package uk.ac.sanger.artemis.util;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 *
 * @author Eric Rasche <rasche.eric@yandex.ru>
 */
public class DatabaseLocationParser {

    private String host;
    private String database;
    private int port = 0;
    private String db_engine = "postgresql";
    private String protocol = "jdbc";
    private static org.apache.log4j.Logger logger4j =
            org.apache.log4j.Logger.getLogger(DatabaseLocationParser.class);
    /**
     * Desire use of Protocol in the final URL. i.e., "jdbc:"
     */
    public static final int PROTOCOL = 1;
    /**
     * Desire use of scheme in final url. i.e., "postgres://"
     */
    public static final int SCHEME = 2;
    /**
     * Desire use of database name in final url.
     *
     */
    public static final int DATABASE_NAME = 4;
    /**
     * Desire listing of query parameters in final url. e.g.,
     * "user=name&ssl=true"
     */
    public static final int QUERY_PARAMS = 8;
    private Map<String, String> params = new HashMap<String, String>();

    /**
     * Empty initializer
     */
    public DatabaseLocationParser() {
    }

    /**
     * Create a new DLP object from a given URL
     *
     * @param url
     */
    public DatabaseLocationParser(String url) {
        setFromURLString(url);
    }

    /**
     * Set the URL internally and parse out important portions.
     *
     * @param url
     */
    private void setFromURLString(String url) {
        logger4j.debug("DLP was called with a URL of [" + url + "]");
        try {
            //"jdbc:postgres://localhost:5432/drupal6?user=drupal6&ssl=true";

            //If it's prefixed, remove that so URI parsing is correct
            if (url.startsWith("jdbc:")) {
                url = url.substring(5);
            }
            if (!url.startsWith(db_engine + "://")) {
                url = db_engine + "://" + url;
            }


            URI db_loc = new URI(url);

            logger4j.debug("URI " + db_loc.toString());
            logger4j.debug("Host: " + db_loc.getHost());
            logger4j.debug("Port: " + db_loc.getPort());
            logger4j.debug("Engine: " + db_loc.getScheme());
            logger4j.debug("DB: " + db_loc.getPath());

            host = db_loc.getHost();

            database = db_loc.getPath().substring(1);

            port = db_loc.getPort();

            db_engine = db_loc.getScheme();

            //Split on '&' and parse each subunit
            String[] query_params = db_loc.getQuery().split("&");
            for (int i = 0; i < query_params.length; i++) {
                //Split based on the equals sign
                logger4j.debug("Given a parameter:" + query_params[i]);
                String[] parts = query_params[i].split("=");


                //This will fail for input like user=chad=o&ssl=true
                //As we'll only grab user=chad
                // Then again, who has an equals sign in their username
                if (parts.length > 1) {
                    params.put(parts[0], parts[1]);
                    logger4j.debug("[" + parts[0] + "," + parts[1] + "]");

                } else {
                    // This might fail strangely, but then again they're providing funky URLs
                    params.put("user", parts[0]);
                    logger4j.debug("[user," + parts[0] + "]");
                }
            }
        } catch (URISyntaxException ex) {
            logger4j.warn("Error parsing URL [" + url + "]" + ex);
        }
        logger4j.debug("This has a complete_url of [" + getCompleteURL() + "]");
    }

    /**
     * Returns the complete URL
     *
     * @return complete url
     */
    public String getCompleteURL() {
        return getURLWithFixes(PROTOCOL | SCHEME
                | DATABASE_NAME | QUERY_PARAMS);
    }

    /**
     * Returns the URL as required for a DriverManager.getConnection object
     *
     * @return url as required for a SQL connection
     */
    public String getConnectionString() {
        return getURLWithFixes(PROTOCOL | SCHEME
                | DATABASE_NAME | QUERY_PARAMS);
    }

    /**
     * Returns the unprefixed URL, for classes that automatically prepend
     * 'jdbc:postgres://'.
     *
     * This is important for uk/ac/sanger/artemis/chado/DbSqlConfig.java
     *
     * @return unprefixed url with database and query parameters appended.
     */
    public String getUnprefixedURL() {
        return getURLWithFixes(DATABASE_NAME | QUERY_PARAMS);
    }

    /**
     * Returns a URL with a selection of modifications.
     *
     * Using a binary OR, one can select which modifications should be applied
     * to create the final URL. These modifications consist of:
     *
     * - PROTOCOL: "jdbc:" - SCHEME: "postgresql://" - DATABASE_NAME: The
     * supplied database name - QUERY_PARAMS: The supplied query parameters
     * (user, ssl, etc)
     *
     * @param modifications, a binary OR'd selection of PROTOCOL, SCHEME,
     * DATABASE_NAME, QUERY_PARAMS, all of which are available as public final
     * static integers from this class
     * @return String version of a URL, modified according to the rules
     * supplied.
     */
    public String getURLWithFixes(int modifications) {
        try {
            String scheme = new String(db_engine);
            String userInfo = null;
            int db_port = new Integer(port);
            String db_name = "/" + database;
            String query_params = "";
            String fragment = null;

            String result = "";
            if ((modifications & PROTOCOL) == PROTOCOL) {
                result += protocol + ":";
            }
            if ((modifications & SCHEME) != SCHEME) {
                scheme = null;
            }
            if ((modifications & DATABASE_NAME) != DATABASE_NAME) {
                db_name = null;
            }

            // Query Parameters
            // "user=chado_user&ssl=true"
            if ((modifications & QUERY_PARAMS) == QUERY_PARAMS) {
                if (params.size() > 0) {
                    /**
                     * Handling of the query parameters. There are other ways to
                     * do this, but it's probably never going to be more than a
                     * "user" and an "ssl" parameter, which means that in the
                     * grand scheme of things, joining a couple strings together
                     * isn't a big issue
                     */
                    Set<String> keys = params.keySet();
                    java.util.Iterator<String> it = keys.iterator();
                    while (it.hasNext()) {
                        String key = it.next();
                        query_params += key + "=" + params.get(key) + "&";
                    }
                    query_params = query_params.substring(0, query_params.length() - 1);
                }
            }
            URI uri_result = new URI(scheme, userInfo, host, db_port, db_name, query_params, fragment);
            logger4j.debug("Pre-final URL: " + uri_result.toString());

            result = result + uri_result.toString();
            //Bugfix. Even if SCHEME is null, // is still prepended. so we remove
            if ((modifications & SCHEME) != SCHEME && (modifications & PROTOCOL) != PROTOCOL) {
                result = result.substring(2);
            }
            return result;
        } catch (URISyntaxException ex) {
            logger4j.error("Could not construct URL. This will likely cause an "
                    + "SQL connection failure. ");
        }
        return null;
    }

    /**
     * Checks whether SSL is enabled
     *
     * @return true if SSL is enabled
     */
    public boolean isSSLEnabled() {
        if (params.containsKey("ssl")) {
            return params.get("ssl").equals("true");
        } else {
            return false;
        }
    }

    /**
     * Returns the hostname
     *
     * @return the hostname
     */
    public String getHost() {
        return host;
    }

    /**
     * Returns the database
     *
     * @return the database name
     */
    public String getDatabase() {
        return database;
    }

    /**
     * Returns the port number
     *
     * @return port number
     */
    public int getPort() {
        return port;
    }

    /**
     * Returns the username of the user connecting to the database
     *
     * @return username
     */
    public String getUsername() {
        if (params.containsKey("user")) {
            return params.get("user");
        } else {
            return "chado"; //That's the default u/n afaik
        }
    }

    /**
     * Returns the database engine
     *
     * @return database engine (usu. postgresql)
     */
    public String getDBEngine() {
        return db_engine;
    }

    /**
     * Sets the hostname
     *
     * @param hostname
     */
    public void setHost(String hostname) {
        host = hostname.trim();
    }

    /**
     * Sets the database name
     *
     * @param new_db_name
     */
    public void setDatabase(String new_db_name) {
        database = new_db_name.trim();
    }

    /**
     * Sets the port number
     *
     * @param new_port_number
     */
    public void setPort(String new_port_number) {
        port = Integer.parseInt(new_port_number.trim());
    }

    /**
     * Sets the port number
     *
     * @param new_port_number
     */
    public void setPort(int new_port_number) {
        port = new_port_number;
    }

    /**
     * Sets the username to connect with
     *
     * @param new_username
     */
    public void setUsername(String new_username) {
        params.put("user", new_username.trim());
    }

    /**
     * Enables or disables SSL in connection URL
     *
     * @param is_enabled
     */
    public void setSSL(boolean is_enabled) {
        if (is_enabled) {
            params.put("ssl", "true");
        } else {
            params.remove("ssl");
        }
    }

    /**
     * Sets the Database Engine
     *
     * @param new_db_engine
     */
    public void setDBEngine(String new_db_engine) {
        db_engine = new_db_engine.trim();
    }
}
