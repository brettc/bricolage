import logging
log = logging.getLogger("experimentdb")

# SQL stuff
from sqlalchemy import (create_engine, MetaData, Table, Column, Integer,
    String, Float, ForeignKey)
from sqlalchemy.sql import and_, or_, not_, delete
from sqlalchemy.orm import mapper, relationship, Session
from sqlalchemy.ext.declarative import declarative_base
SQLBase = declarative_base()


class TreatmentRecord(SQLBase):
    __tablename__ = 'treatment'
    treatment_id = Column(Integer, primary_key=True)
    name = Column(String(30))
    replicates = relationship("ReplicateRecord", cascade="all, delete-orphan",
                               backref='treatment')

    def __init__(self, t):
        self.treatment_id = t.seq
        self.name = t.name

    def __str__(self):
        return "Treatment:{0.treatment_id:03d}, {0.name:}".format(self)

class ReplicateRecord(SQLBase):
    __tablename__ = 'replicate'
    treatment_id = Column(Integer, ForeignKey('treatment.treatment_id'), 
                          primary_key=True)
    replicate_id = Column(Integer, primary_key=True)
    seed = Column(Integer)
    generations = Column(Integer)

    def __init__(self, r):
        self.treatment_id = r.treatment.seq
        self.replicate_id = r.seq
        self.generations = r.generations
        self.seed = r.seed

    def __str__(self):
        return "Replicate:{0.treatment_id:02d}:{0.replicate_id:02d}, "\
            "S{0.seed:010d}, G{0.generations:010d}".format(self)


class Database(object):
    def __init__(self, folder):
        self.path = (folder / 'stats').with_suffix('.sqlite')
        self.engine = None
        self.session = None

    def create(self, echo=False):
        log.info("Initialising database at {}".format(str(self.path)))
        self.engine = create_engine('sqlite:///{}'.format(str(self.path)), 
                                    echo=echo)
        SQLBase.metadata.create_all(self.engine)
        self.session = Session(self.engine)

    def save_frame(self, table, replicate, frame):
        if self.engine.dialect.has_table(self.engine.connect(), table):
            # Delete any existing 
            self.engine.execute("delete from {}"
                                " where treatment_id = {}"
                                " and replicate_id = {}".format(
                table, replicate.treatment.seq, replicate.seq))

        # Add the treatment / replicate columns
        frame['replicate_id'] = replicate.seq
        frame['treatment_id'] = replicate.treatment.seq

        # Set an index, and add it -- 
        # TODO: Still missing more indexing info
        frame.set_index(['treatment_id', 'replicate_id'], inplace=True)
        frame.to_sql(table, self.engine, if_exists='append', index_label=
                     ['treatment_id', 'replicate_id'])
        

