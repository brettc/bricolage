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

    def create(self, args):
        log.info("Initialising database at {}".format(str(self.path)))
        self.engine = create_engine('sqlite:///{}'.format(str(self.path)), 
                                    echo=args.verbose)
        SQLBase.metadata.create_all(self.engine)
        self.session = Session(self.engine)

