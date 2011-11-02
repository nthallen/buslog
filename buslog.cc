#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "curllog/curl_obj.h"
#include "nortlib.h"

static int done = 0;

const char *find_last_time(xmlNodePtr tree) {
  const char *lastTime = NULL;
  for ( ; tree != NULL; tree = tree->next ) {
    if ( tree->name != NULL && strcmp((const char *)tree->name, "lastTime") == 0 ) {
      lastTime = (const char *)xmlGetProp(tree, (const xmlChar *)"time" );
      return lastTime;
    }
    lastTime = find_last_time( tree->children );
    if ( lastTime != NULL ) return lastTime;
  }
  return NULL;
}

static const char *get_prop( xmlNodePtr tree, const char *prop ) {
  const char *pval = (const char *)xmlGetProp(tree, (const xmlChar *)prop);
  if (pval == NULL) pval = "n/a";
  return pval;
}

void write_log_int( FILE *fp, const char *lastTime, xmlNodePtr tree ) {
  for ( ; tree != NULL; tree = tree->next ) {
    if ( tree->name != NULL && strcmp((const char *)tree->name, "vehicle") == 0) {
      const char *id = get_prop( tree, "id" );
      const char *dirTag = get_prop(tree, "dirTag");
      const char *lat = get_prop(tree, "lat");
      const char *lon = get_prop(tree, "lon");
      const char *secsSinceReport = get_prop(tree, "secsSinceReport");
      const char *heading = get_prop(tree, "heading");
      fprintf( fp, "%s, %s, \"%s\", \"%s\", %s, %s, %s\n",
        lastTime, secsSinceReport, id, dirTag, lat, lon, heading );
    }
    write_log_int( fp, lastTime, tree->children );
  }
}

void write_log( const char *logfile, const char *lastTime, xmlNodePtr tree ) {
  FILE *fp = fopen( logfile, "a" );
  if ( fp == NULL )
    nl_error(4, "Could not append to log file %s", logfile );
  write_log_int( fp, lastTime, tree );
  fclose(fp);
}

int main( int argc, char **argv ) {
  const char *lastTime = "0";
  const char *route = "77";
  const char *logfile = "locations.log";
  curl_obj co;
  co.set_log_level( CT_LOG_BODIES );
  while (!done) {
    xmlNodePtr tree;
    char urlbuf[256];
    snprintf( urlbuf, 256,
      "http://webservices.nextbus.com/service/publicXMLFeed?command=vehicleLocations&a=mbta&r=%s&t=%s",
      route, lastTime );
    co.transaction_start("Get bus data");
    co.set_url(urlbuf);
    co.perform("Get Locations");
    co.transaction_end();
    tree = co.get_parse_tree();
    lastTime = find_last_time(tree);
    if (lastTime != NULL) {
      write_log( logfile, lastTime, tree );
    } else lastTime = "0";
    co.write_transaction_log();
    sleep(30);
  }
  return 0;
}
