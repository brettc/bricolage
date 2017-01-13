import logging

from sqlalchemy import (create_engine, Column, Integer,
                        String, Float, ForeignKey)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, Session

log = logging.getLogger("experimentdb")

SQLBase = declarative_base()


class TreatmentRecord(SQLBase):
    __tablename__ = 'treatment'
    treatment_id = Column(Integer, primary_key=True)
    name = Column(String(30))
    # seed = Column(Integer)
    replicates = relationship("ReplicateRecord", cascade="all, delete-orphan",
                              backref='treatment')

    def __init__(self, t):
        self.treatment_id = t.seq
        self.name = t.name
        # self.seed = t.seed

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
        return "Replicate:{0.treatment_id:02d}:{0.replicate_id:02d}, " \
               "S{0.seed:010d}, G{0.generations:010d}".format(self)


class StatsRecord(SQLBase):
    __tablename__ = 'stats'
    treatment_id = Column(Integer,
                          ForeignKey('treatment.treatment_id'),
                          primary_key=True)
    replicate_id = Column(Integer,
                          ForeignKey('replicate.replicate_id'),
                          primary_key=True)
    generation = Column(Integer, primary_key=True)
    kind = Column(String(10), primary_key=True)
    tag = Column(String(10), index=True)
    value = Column(Float())

    def __init__(self, rep, generation_number, kind, val, tag):
        self.treatment_id = rep.treatment.seq
        self.replicate_id = rep.seq
        self.generation = generation_number
        self.kind = kind
        self.value = val
        self.tag = tag


class StatsGroupRecord(SQLBase):
    __tablename__ = 'stats_gen'
    treatment_id = Column(Integer,
                          ForeignKey('treatment.treatment_id'),
                          primary_key=True)
    replicate_id = Column(Integer,
                          ForeignKey('replicate.replicate_id'),
                          primary_key=True)
    tag = Column(String(10), primary_key=True)
    generation = Column(Integer, primary_key=True)

    def __init__(self, rep, generation_number, tag):
        self.treatment_id = rep.treatment.seq
        self.replicate_id = rep.seq
        self.generation = generation_number
        self.tag = tag


class StatsReplicateRecord(SQLBase):
    __tablename__ = 'stats_rep'
    treatment_id = Column(Integer,
                          ForeignKey('treatment.treatment_id'),
                          primary_key=True)
    replicate_id = Column(Integer,
                          ForeignKey('replicate.replicate_id'),
                          primary_key=True)
    kind = Column(String(10), primary_key=True)
    value = Column(Float())

    def __init__(self, rep, kind, val):
        self.treatment_id = rep.treatment.seq
        self.replicate_id = rep.seq
        self.kind = kind
        self.value = val


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

    def drop_table(self, table):
        self.engine.execute("drop table if exists {}".format(table))

    def save_frame(self, table, replicate, frame):
        if self.engine.dialect.has_table(self.engine.connect(), table):
            # Delete any existing
            self.engine.execute(
                "delete from {}"
                " where treatment_id = {}"
                " and replicate_id = {}".format(
                    table, replicate.treatment.seq, replicate.seq))

        # Add the treatment / replicate columns
        frame['replicate_id'] = replicate.seq
        frame['treatment_id'] = replicate.treatment.seq

        # Set an index, and add it --
        # TODO: Still missing more indexing info
        frame.set_index(['treatment_id', 'replicate_id'], inplace=True)
        frame.to_sql(table,
                     self.engine,
                     if_exists='append',
                     index_label=['treatment_id', 'replicate_id'])

    def remove(self, rep, kind=None):
        template = "delete from stats where "
        template += " treatment_id = {} and replicate_id = {}"
        args = [rep.treatment.seq, rep.seq]
        if kind is not None:
            template += " and kind = '{}'"
            args.append(kind)
        self.engine.execute(template.format(*args))
