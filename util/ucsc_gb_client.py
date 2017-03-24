from sqlalchemy import create_engine
from sqlalchemy.engine.url import URL
from sqlalchemy.pool import NullPool
import pandas as pd
import re

chia_pet_tables = dict(
    Hct116Pol2Rep1="wgEncodeGisChiaPetHct116Pol2InteractionsRep1",
    Helas3Pol2Rep1="wgEncodeGisChiaPetHelas3Pol2InteractionsRep1",
    K562CtcfRep1="wgEncodeGisChiaPetK562CtcfInteractionsRep1",
    K562Pol2Rep1="wgEncodeGisChiaPetK562Pol2InteractionsRep1",
    K562Pol2Rep2="wgEncodeGisChiaPetK562Pol2InteractionsRep2",
    Mcf7CtcfRep1="wgEncodeGisChiaPetMcf7CtcfInteractionsRep1",
    Mcf7CtcfRep2="wgEncodeGisChiaPetMcf7CtcfInteractionsRep2",
    Mcf7EraaRep1="wgEncodeGisChiaPetMcf7EraaInteractionsRep1",
    Mcf7EraaRep2="wgEncodeGisChiaPetMcf7EraaInteractionsRep2",
    Mcf7EraaRep3="wgEncodeGisChiaPetMcf7EraaInteractionsRep3",
    Mcf7Pol2Rep1="wgEncodeGisChiaPetMcf7Pol2InteractionsRep1",
    Mcf7Pol2Rep2="wgEncodeGisChiaPetMcf7Pol2InteractionsRep2",
    Mcf7Pol2Rep3="wgEncodeGisChiaPetMcf7Pol2InteractionsRep3",
    Mcf7Pol2Rep4="wgEncodeGisChiaPetMcf7Pol2InteractionsRep4",
    Nb4Pol2Rep1="wgEncodeGisChiaPetNb4Pol2InteractionsRep1",
)


class UcscGbClient:
    __db_url = dict(
        drivername='mysql+pymysql',
        host='genome-mysql.cse.ucsc.edu',
        port='3306',
        username='genome',
        password='',
        database='hg19',
        query={'charset': 'utf8'}
    )

    def __init__(self):
        # db = create_engine('mysql://bud:earth@localhost:3306/hg19') # require module `MySQLdb`
        #   default dialect is 'mysql+mysql-python'
        #   `MySQLdb` is a fork of MySQL-python with added support for Python 3
        #   See http://docs.sqlalchemy.org/en/latest/core/engines.html#mysql

        # db = create_engine('mysql+pymysql://bud:earth@localhost:3306/hg19') # require module `PyMySQL`

        # For `poolclass`, see http://stackoverflow.com/a/8705750

        self.db = create_engine(URL(**UcscGbClient.__db_url), poolclass=NullPool)
        self.conn = self.db.connect()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
        self.db.dispose()

    def query_chia_pet_cluster(self, table_key, scope_type, chrom, chrom_start, chrom_end):
        if table_key not in chia_pet_tables:
            raise ValueError("No such table key: '{}'".format(table_key))

        if scope_type == 'within':
            scope_condition = "chrom='{q_chrom}' and chromStart>={q_chrom_start} and chromEnd<={q_chrom_end}".\
                format(q_chrom=chrom, q_chrom_start=chrom_start, q_chrom_end=chrom_end)
        elif scope_type == 'overlap':
            scope_condition = "chrom='{q_chrom}' and chromStart<{q_chrom_end} and chromEnd>{q_chrom_start}".\
                format(q_chrom=chrom, q_chrom_start=chrom_start, q_chrom_end=chrom_end)
        else:
            raise ValueError("Scope type '{}' is not supported".format(scope_type))

        query = '''
                SELECT chrom, chromStart, chromEnd, name
                FROM {tableName}
                WHERE {condition}
                '''.format(tableName=chia_pet_tables[table_key], condition=scope_condition)

        rows = self.conn.execute(query)

        data = rows.fetchall()
        df = pd.DataFrame(data)

        if df.empty:
            return df

        df.columns = ['chrom', 'chromStart', 'chromEnd', 'blocks']
        id_col = df.apply(lambda x: "{}_{}_{}_{}".format(table_key, x['chrom'], x['chromStart'], x['chromEnd']), axis=1)

        df_block = pd.DataFrame(list(map(lambda cell: re.split(":|\.\.|-|,", cell)[:-1], df.loc[:, 'blocks'])))
        df_block.columns = ['b1_chrom', 'b1_chromStart', 'b1_chromEnd', 'b2_chrom', 'b2_chromStart', 'b2_chromEnd']
        for col in ['b1_chromStart', 'b1_chromEnd', 'b2_chromStart', 'b2_chromEnd']:
            df_block.loc[:, col] = pd.to_numeric(df_block[col])

        df_block = df_block.assign(id=id_col)

        return df_block
