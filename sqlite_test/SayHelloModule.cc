// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include <arcane/core/Directory.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
  Directory d(subDomain()->exportDirectory(), "test.db");

  Integer rc = sqlite3_open(d.path().localstr(), &m_db);
  if(rc){
    sqlite3_close(m_db);
    ARCANE_FATAL("Can't open database : {0}", sqlite3_errmsg(m_db));
  }
}

void SayHelloModule::
createTable()
{
  info() << "Module SayHello INIT";
  char *zErrMsg = nullptr;
  Integer rc = sqlite3_exec(m_db, "CREATE TABLE IF NOT EXISTS utilisateur(id INT, truc TEXT)", nullptr, nullptr, &zErrMsg);

  if( rc != SQLITE_OK ){
    error() << "SQL error : " << zErrMsg;
    sqlite3_free(zErrMsg);
  }
}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";
  char *zErrMsg = nullptr;
  Integer rc = sqlite3_exec(m_db, "INSERT INTO utilisateur (id, truc) VALUES (123, \"abcd\")", nullptr, nullptr, &zErrMsg);

  if( rc != SQLITE_OK ){
    error() << "SQL error : " << zErrMsg;
    sqlite3_free(zErrMsg);
  }
}

void SayHelloModule::
printResults()
{
  info() << "Module SayHello END";

  sqlite3_stmt *stmt;

  sqlite3_prepare_v2(m_db, "SELECT * FROM utilisateur", -1, &stmt, nullptr);

  Integer step = sqlite3_step(stmt);

  while (step != SQLITE_DONE && step != SQLITE_MISUSE) {
    for (Integer i = 0; i < sqlite3_column_count(stmt); ++i) {
      info() << "Column : " << i;
      switch (sqlite3_column_type(stmt, i))
      {
      case (SQLITE_INTEGER):
        info() << "SQLITE_INTEGER : " << sqlite3_column_int(stmt, i);
        break;
      case (SQLITE_FLOAT):
        info() << "SQLITE_FLOAT : " << sqlite3_column_double(stmt, i);
        break;
      case (SQLITE3_TEXT): {
        Integer size = sqlite3_column_bytes(stmt, i);
        info() << "size of text : " << size;
        String str(reinterpret_cast<const char*>(sqlite3_column_text(stmt, i)));
        String str1(str.clone());
        info() << "size of text 1 : " << str1.length();
        info() << "SQLITE3_TEXT : " << str1;
        break;
      }
    case (SQLITE_BLOB): {
        info() << "SQLITE_BLOB (Binary Large Object)";
        Integer size = sqlite3_column_bytes(stmt, i);
        Byte* blob = (Byte*)sqlite3_column_blob(stmt, i);
        break;
      }
      case (SQLITE_NULL):
        info() << "SQLITE_NULL";
        break;
      default:
        ARCANE_FATAL("Unknown type");
      }
    }
    step = sqlite3_step(stmt);
  }
  sqlite3_finalize(stmt);
}

void SayHelloModule::
endModule()
{
  info() << "Module SayHello END";
  sqlite3_close(m_db);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SAYHELLO(SayHelloModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
