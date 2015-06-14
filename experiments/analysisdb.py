
# SQL stuff
from sqlalchemy import (create_engine, MetaData, Table, Column, Integer,
    String, Float, ForeignKey)
from sqlalchemy.sql import and_, or_, not_, delete
from sqlalchemy.orm import mapper, relationship, Session
from sqlalchemy.ext.declarative import declarative_base
SQLBase = declarative_base()

# --------------------------------------------
# Define the reporting objects
# _stat_types = 'mean', 'count'
# _classes = Sender, Receiver, Loner, Modular
# _names = [cls.__name__ for cls in _classes]
#

class Replicate(SQLBase):
    __tablename__ = 'replicate'
    replicate_id = Column(Integer, primary_key=True)
    treatment = Column(String(40))
    rep_num = Column(Integer)
    stats = relationship("Stats", cascade="all, delete-orphan",
                               backref='treatment')

    def __init__(self, treatment, rep):
        self.treatment = treatment
        self.rep_num = rep

    def __repr__(self):
        return "<Replicate: {0.treatment} / {0.rep_num}>".format(self)

class Stats(SQLBase):
    __tablename__ = 'stats'

    replicate_id = Column(Integer, ForeignKey('replicate.replicate_id'), 
                          primary_key=True)
    generation = Column(Integer, primary_key=True)
    kind = Column(String(10), primary_key=True)
    value = Column(Float)

    def __init__(self, generation, kind, value):
        self.generation = generation
        self.kind = kind
        self.value = value


class AnalysisDB(object):
    def __init__(self, path, verbose=False):
        self.engine = create_engine('sqlite:///{}'.format(self.path), echo=verbose)
        SQLBase.metadata.create_all(self.engine)
        self.session = Session(self.engine)
