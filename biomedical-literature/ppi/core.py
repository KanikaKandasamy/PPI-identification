# -*- encoding: latin1 -*-

import requests
from xml.dom import minidom
import json
import re
import time
from pycorenlp import StanfordCoreNLP
import itertools
from bisect import bisect_left
import math
import sys
import pkg_resources
from scipy import sparse
import logging
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

try:
    import cPickle as pickle
    from HTMLParser import HTMLParser
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    import html
    import _pickle as pickle
    from importlib import reload


NLP = StanfordCoreNLP('http://localhost:9000')

def pmid_2_pmc(identifiers):
    pmcids = set()
    maxidents = 200

    for subset in [identifiers[x:x+maxidents] for x in range(0, len(identifiers),maxidents)]:
        params = {
            'ids': ",".join(subset),
            'format': 'json'
        }
        req = requests.get("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/", params=params)
        if req.status_code == 200:
            response = json.loads(req.content.decode('latin1'))
            for record in response['records']:
                if 'status' in record:
                    continue
                pmcids.add(record['pmcid'][3:])
        else:
            raise PubMedQueryError("Can't convert identifiers through Pubmed idconv tool.")
        time.sleep(3)
    return list(pmcids)

def take_closest(mylist, mynumber):
    pos = bisect_left(mylist, mynumber)
    if pos == 0:
        return mylist[0]
    if pos == len(mylist):
        return mylist[-1]
    before = mylist[pos - 1]
    after = mylist[pos]
    if after - mynumber < mynumber - before:
        return after
    else:
        return before

def take_farthest(mylist, mynumber):
    distance_first = mylist[0]  - mynumber
    distance_last  = mylist[-1] - mynumber
    if math.fabs(distance_first) > math.fabs(distance_last):
        return mylist[0]
    else:
        return mylist[-1]

def minidom_to_text(minidom):
    return " ".join(t.nodeValue for t in minidom.childNodes if t.nodeType == t.TEXT_NODE)


class PMQuery(object):
    def __init__(self, ids, database="PMC"):
        self.ids = ids
        self.database = database
        self.articles = list()
        self.found    = set()
        self.notfound = set()

    def __get_pmc(self, req):
        if req.status_code == 200:
            article_text = minidom.parseString(req.content)
            articles = article_text.getElementsByTagName('article')
            for article in articles:
                pmid    = article.getElementsByTagName('article-id')[0].firstChild.nodeValue
                journal = article.getElementsByTagName('journal-id')[0].firstChild.nodeValue
                year = article.getElementsByTagName('year')[0].firstChild.nodeValue
                try:
                    pmcid = article.getElementsByTagName('article-id')[1].firstChild.nodeValue
                except:
                    continue
                body =  article.getElementsByTagName('body')
                if len(body) == 0:
                    continue
                self.found.add(pmid)
                paragraphs = body[0].getElementsByTagName('p')
                fulltext = list()
                for par in paragraphs:
                    fulltext.append(minidom_to_text(par))
                self.articles.append(Article(pmid=pmid, pmcid=pmcid, journal=journal, year=year, fulltext="\n".join(fulltext)))
            self.notfound = set(self.ids).difference(self.found)
        else:
            PubMedQueryError("Can't connect to PMC...")

    def __get_pubmed(self, req):
        if req.status_code == 200:
            article_text = minidom.parseString(req.content)
            articles = article_text.getElementsByTagName('PubmedArticle')
            for article in articles:
                pmid = article.getElementsByTagName('PMID')[0]
                pmid_text = minidom_to_text(pmid)
                journal = article.getElementsByTagName('Journal')[0].getElementsByTagName('Title')[0]
                journal_text = minidom_to_text(journal)
                year = article.getElementsByTagName('Year')[0].firstChild.nodeValue
                abstracts = article.getElementsByTagName('AbstractText')
                abstract_text = list()
                for abst in abstracts:
                    abstract_text.append(minidom_to_text(abst))
                abstract_text = "\n".join(abstract_text)
                if not abstract_text.strip():
                    continue
                self.found.add(pmid_text)
                self.articles.append(Article(pmid=pmid_text, journal=journal_text, year=year, abstract=abstract_text))
            self.notfound = set(self.ids).difference(self.found)
        else:
            PubMedQueryError("Can't connect to PubMed...")

    def get_articles(self):
        maxidents = 200 

        for subset in [self.ids[x:x+maxidents] for x in range(0, len(self.ids), maxidents)]:
            if self.database == "PMC":

                params = {
                    'id': ",".join(pmid_2_pmc(subset)),
                    'db': 'pmc',
                }
                req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
                self.__get_pmc(req)
            elif self.database == "PUBMED":
                params = {
                    'id':      ",".join(subset),
                    'db':      'pubmed',
                    'retmode': 'xml'
                }
                req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
                self.__get_pubmed(req)
            else:
                logging.error('%s: Incorrect database. Choose "PMC" or "PUBMED"', self.database)

    def __iter__(self):
        return iter(self.articles)

    def __getitem__(self, index):
        return self.articles[index]

# ----------------------------------------------
class Article(object):
    
    def __init__(self, pmid, pmcid=None, journal=None, year=None, fulltext=None, abstract=None):
        self.pmid       = pmid
        self.pmcid      = pmcid
        self.journal    = journal
        self.year       = year
        self.abstract   = abstract
        self.fulltext   = fulltext
        self.sentences  = list()

    def extract_interactions(self, source="fulltext", only_dict=False):
        self.extract_sentences(source=source)
        for sentence in self.sentences:
            sentence.extract_interactions(only_dict)
        

    @property
    def predictions(self):
        predictions = []
        for sentence in self.sentences:
            predictions.extend(sentence.get_interactions())
        return predictions
        
    def as_html(self):
        pass

    def extract_sentences(self, mode="split", source="fulltext"):
        text = ""
        if source == "fulltext":
            text = str(self.fulltext)
        else:
            text = str(self.abstract)

        if mode == "no-split":
            self.sentences.append(Sentence(originaltext=text))
        else:
            caps = "([A-Z])"
            prefixes = "(Mr|Fig|fig|St|Mrs|Ms|Dr)[.]"
            digits = "([0-9])"
            fig_letters = "([A-Ka-k])"
            suffixes = "(Inc|Ltd|Jr|Sr|Co)"
            starters = r"(Mr|Mrs|Ms|Dr|He\s|She\s|It\s|They\s|Their\s|Our\s|We\s|But\s|However\s|That\s|This\s|Wherever)"
            acronyms = "([A-Z][.][A-Z][.](?:[A-Z][.])?)"
            websites = "[.](com|net|org|io|gov)"
            species = r"([A-Z])[.] ?([a-z]+)"
            text = " " + text + "  "
            text = text.replace("\n"," ")
            text = re.sub(prefixes,"\\1<prd>",text)
            text = re.sub(websites,"<prd>\\1",text)
            if "Ph.D" in text:
                text = text.replace("Ph.D.","Ph<prd>D<prd>")
            text = re.sub(r"\s" + caps + "[.] "," \\1<prd> ",text)
            text = re.sub(acronyms+" "+starters,"\\1<stop> \\2",text)
            text = re.sub(caps + "[.]" + caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>\\3<prd>",text)
            text = re.sub(caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>",text)
            text = re.sub(" "+suffixes+"[.] "+starters," \\1<stop> \\2",text)
            text = re.sub(" "+suffixes+"[.]"," \\1<prd>",text)
            text = re.sub(" " + caps + "[.]"," \\1<prd>",text)
            text = re.sub(digits + caps + "[.]"," \\1<prd>",text)
            text = re.sub(digits + "[.]" + digits,"\\1<prd>\\2",text)
            text = re.sub(digits + "[.]" + fig_letters,"\\1<prd>\\2",text)
            text = re.sub(species, "\\1<prd> \\2", text)
            if "”"    in text:
                text = text.replace(".”","”.")
            if "\""   in text:
                text = text.replace(".\"","\".")
            if "!"    in text:
                text = text.replace("!\"","\"!")
            if "?"    in text:
                text = text.replace("?\"","\"?")
            if "e.g." in text:
                text = text.replace("e.g.","e<prd>g<prd>")
            if "i.e." in text:
                text = text.replace("i.e.","i<prd>e<prd>")
            text = text.replace(".",".<stop>")
            text = text.replace("?","?<stop>")
            text = text.replace("!","!<stop>")
            text = text.replace("<prd>",".")
            sentences = text.split("<stop>")
            #sentences = sentences[:-1]
            sentences = [ s.strip() for s in sentences ]
            try:
                html = HTMLParser()
            except:
                import html
            for sentence in sentences:
                sentence = str(html.unescape(sentence))
                if not sentence.strip() or not isinstance(sentence, str):
                    continue
                self.sentences.append(Sentence(originaltext=sentence))

    def count_genes(self):
        pass

    def make_wordcloud(self):
        pass

    def __str__(self):
        return "Article with PubMED id:%s" % (self.pmid)

# ----------------------------------------------
class Protein(object):
    GENEDICT = dict()
    GENEDICTFILE = pkg_resources.resource_filename('ppi', 'data/HGNC_gene_dictionary.txt')
    try:
        with open(GENEDICTFILE, 'r') as f:
            for line in f:
                line = line.strip()
                cols = line.split("\t")
                GENEDICT[cols[0]] = cols[0]
                for alias in cols[1:]:
                    GENEDICT[alias.upper()] = cols[0]
    except Exception:
        raise FileNotFoundError("Can't read gene dictionary file %s\n" % GENEDICTFILE)


    def __init__(self, symbol, positions, sentence):
        self.symbol = symbol
        self.positions = positions
        self.sentence = sentence
        self.synonym = list()
        self.count = len(positions)

    def is_in_dict(self):
        disambiguated = self.symbol.upper()
        disambiguated = disambiguated.replace("'", "")
        disambiguated = disambiguated.replace('"', '')
        if disambiguated in Protein.GENEDICT:
            return True
        else:
            return False

    def disambiguate(self):
        disambiguated = self.symbol.upper()
        disambiguated = disambiguated.replace("'", "")
        disambiguated = disambiguated.replace('"', '')
        if disambiguated in Protein.GENEDICT:
            return Protein.GENEDICT[disambiguated]
        else:
            return disambiguated

    def __str__(self):
        return "%s found in positions %s" % (self.symbol, ":".join([ str(idx) for idx in self.positions ]))


# ----------------------------------------------
class Sentence(object):
    def __init__(self, originaltext):
        self.originaltext = originaltext
        self.tokens       = list()
        self.tree         = list()
        self.candidates   = list()
        self.proteins     = list()

    def annotate(self):
        if not self.originaltext.strip():
            self.tokens = ""
        annotated = json.loads(NLP.annotate(self.originaltext))
        if annotated['sentences']:
            self.tokens = annotated['sentences'][0]['tokens']

    def get_candidates(self, only_dict=False):
        if not self.tokens:
            self.annotate()
        if self.candidates:
            return
            
        prot_counter = 0
        state        = 0
        prot_list = list()
        for token in self.tokens:
            if token['ner'] == "P":
                state = 1
                if len(prot_list) - 1 < prot_counter :
                    prot_list.append(list())
                prot_list[prot_counter].append(token['index'])
            else:
                if state == 1:
                    prot_counter += 1
                    state = 0
        prots_in_sentence = list()
        for prot_pos in prot_list:
            protein_symbol = " ".join([self.tokens[token_idx - 1]['word'] for token_idx in prot_pos])
            protein = Protein(
                symbol=protein_symbol,
                positions=prot_pos,
                sentence=self
            )
            if only_dict is True:
                if not protein.is_in_dict():
                    continue
            self.proteins.append(protein)
            prots_in_sentence.append(protein)

        for prot in itertools.combinations(prots_in_sentence, r=2):
            self.candidates.append(InteractionCandidate(prot1=prot[0], prot2=prot[1]))

    def to_html(self):
        if not self.tokens:
            self.annotate()
        if not self.candidates:
            self.get_candidates()

        state = 0
        html_list = list()
        for token in self.tokens:
            word = token['word']
            word = re.sub("-LRB-", "(", word)
            word = re.sub("-RRB-", ")", word)
            if state == 0:
                if token['ner'] == "O":
                    if re.match('VB[DGNPZ]?', token['pos']):
                        html_list.append('<span class="verb">%s</span>' % word)
                    else:
                        html_list.append(word)
                else:
                    state = 1
                    html_list.append('<span class="prot">')
                    html_list.append(word)
            else:
                if token['ner'] == "O":
                    html_list.append("</span>")
                    if re.match('VB[DGNPZ]?', token['pos']):
                        html_list.append('<span class="verb">%s</span>' % word)
                    else:
                        html_list.append(word)
                    state = 0
                else:
                    html_list.append(word)
        return " ".join(html_list)

    def extract_interactions(self, only_dict):
        self.annotate()
        self.get_candidates(only_dict)
        for candidate in self.candidates:
            candidate.predict()

    def get_interactions(self):
        interactions = [ candidate for candidate in self.candidates if candidate.label is True ]
        return interactions
        

    def __str__(self):
        return self.originaltext


# ----------------------------------------------
class InteractionCandidate(object):
    verb_scores = dict({
        "acetylate":1, "acylate":1, "amidate":1, "brominate":1, "biotinylate":1,
        "carboxylate":1, "cysteinylate":1, "farnesylate":1, "formylate":1, "hydroxilate":1,
        "hydroxylate":1, "methylate":1, "demethylate":1, "myristoylate":1, "myristylate":1,
        "palmitoylate":1, "palmitylate":1, "phosphorylate":1, "dephosphorylate":1, "pyruvate":1,
        "nitrosylate":1, "sumoylate":1, "ubiquitinylate":1, "ubiquitinate":1, "dissociate":1,
        "assemble":1, "attach":1, "bind":2, "complex":1, "contact":1, "couple":1, "multimerize":1,
        "multimerise":1, "dimerize":1, "dimerise":1, "interact":4, "precipitate":1, "regulate":1,
        "inhibit":1, "activate":3, "downregulate":2, "down-regulate":2, "suppress":2, "upregulate":2,
        "up-regulate":1, "block":1, "inactivate":1, "induce":1, "modify":1, "overexpress":1, "promote":1,
        "stimulate":1, "substitute":1, "catalyze":1, "cleave":1, "conjugate":1, "disassemble":1,
        "discharge":1, "mediate":1, "modulate":1, "repress":1, "transactivate":1
    })

    PRED_FILE = pkg_resources.resource_filename('ppi', 'data/RF_scikit.pkl')
    with open(PRED_FILE, 'rb') as f:
        try:
            predictor = pickle.load(f)
            try:
                predictor = pickle.load(f, encoding='utf-8')
            except Exception as err:
                raise(CannotLoadClassifier("RF_scikit.pkl can't be loaded - %s" % err))



    def __init__(self, prot1, prot2):
        self.prot1 = prot1
        self.prot2 = prot2
        self.between_idxes = (prot1.positions[-1], prot2.positions[0] - 1)
        self.label    = None
        self.votes    = None
        self.feat_cols = list() 
        self.feat_current_col = 0
        self.feat_vals = list() 
        self.features_sparse = None

    def compute_features(self):
        self.__token_distance()
        self.__total_tokens()
        self.__verb_features("between")
        self.__verb_features("all")
        self.__pos_features("between")
        self.__pos_features("all")
        self.__prot_count("between")
        self.__prot_count("all")
        self.__keyword_count("between")

        rows = [0 for val in self.feat_vals]
        self.features_sparse = sparse.coo_matrix((self.feat_vals, (rows, self.feat_cols)), shape=(1,178))

    def __token_distance(self):
        subsentence = self.prot1.sentence.tokens[self.between_idxes[0]:self.between_idxes[1]]
        if len(subsentence) > 0:
            self.feat_cols.append(self.feat_current_col)
            self.feat_vals.append(len(subsentence))
        self.feat_current_col += 1

    def __total_tokens(self):
        if len(self.prot1.sentence.tokens) > 0:
            self.feat_cols.append(self.feat_current_col)
            self.feat_vals.append(len(self.prot1.sentence.tokens))
        self.feat_current_col += 1

    def __prot_count(self, mode="all"):
        prota_count = 0
        protb_count = 0
        tokens = list()
        if mode == "all":
            tokens = self.prot1.sentence.tokens
        else:
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1
            tokens = self.prot1.sentence.tokens[init_coord:final_coord]

        for token in tokens:
            if token['word'] == self.prot1.symbol:
                prota_count += 1
            if token['word'] == self.prot2.symbol:
                protb_count += 1

        self.feat_cols.append(self.feat_current_col)
        self.feat_current_col += 1
        self.feat_cols.append(self.feat_current_col)
        self.feat_vals.extend([prota_count, protb_count])
        self.feat_current_col += 1

    def __get_token_pos(self, mode="all"):
        init_coord  = None
        final_coord = None
        if mode == "all":
            init_coord  = 0
            final_coord = len(self.prot1.sentence.tokens)
        else:
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1

        pos_str = list()
        subsentence = self.prot1.sentence.tokens[init_coord:final_coord]
        for token in subsentence:
            pos_str.append(token['pos'])
        return ",".join(pos_str)

    def __pos_features(self, mode):
        pos_counts = {
            "CC": 0, "LS": 0, "MD": 0,
            "NN": 0, "NNS": 0, "NNP": 0,
            "NNPS": 0, "PDT": 0, "POS": 0,
            "PRP": 0, "PRP$": 0, "RB": 0,
            "RBR": 0, "RBS": 0, "RP": 0,
            "SYM": 0, "TO": 0, "UH": 0,
            "VB": 0, "VBD": 0, "VBG": 0,
            "VBN": 0, "VBP": 0, "VBZ": 0,
            "WDT": 0, "WP": 0, "WP$": 0,
            "WRB": 0, "IN": 0, "DT": 0, ".": 0,
            "JJ" :0, "CD": 0, "": 0, "-LRB-": 0,
            "-RRB-": 0, "JJR": 0, ":": 0, "FW": 0,
            'JJS': 0, "EX": 0, "''": 0, ",":0
        }

        if mode == "all" or mode == "between":
            for pos in self.__get_token_pos(mode=mode).split(","):
                if pos in pos_counts:
                    pos_counts[pos] += 1
        for postag in sorted(pos_counts):
            if pos_counts[postag] > 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(pos_counts[postag])
            self.feat_current_col += 1

    def __verb_distances(self, pidx, vidxes):
        closest_verb_idx  = take_closest(vidxes, pidx)
        farthest_verb_idx = take_farthest(vidxes, pidx)

        closest_distance  = math.fabs(closest_verb_idx -  pidx)
        farthest_distance = math.fabs(farthest_verb_idx -  pidx)
        return (closest_distance, farthest_distance)


    def __verb_features(self, flag):
        numverbs = dict({
            'VB': 0,  'VBD': 0, 'VBG': 0,
            'VBN': 0, 'VBP': 0, 'VBZ': 0
        })
        maxscore   = 0
        totalscore = 0
        verb_idxes = list()
        someverb_flag = False
        tokens_to_work = list()

        if flag == "between":
            tokens_to_work = self.prot1.sentence.tokens[self.between_idxes[0]:self.between_idxes[1]]
        else:
            tokens_to_work = self.prot1.sentence.tokens

        for token in tokens_to_work:
            if re.match('VB[DGNPZ]?', token['pos']):
                someverb_flag = True
                numverbs[token['pos']]  += 1
                verb_idxes.append(token['index'])
                if token['lemma'] in InteractionCandidate.verb_scores:
                    totalscore +=  InteractionCandidate.verb_scores[token['lemma']]
                    if InteractionCandidate.verb_scores[token['lemma']] > maxscore:
                        maxscore = InteractionCandidate.verb_scores[token['lemma']]

        (cl1, far1, cl2, far2) = (0,0,0,0)
        if someverb_flag is True:
            (cl1, far1) = self.__verb_distances(self.prot1.positions[-1], verb_idxes)
            (cl2, far2) = self.__verb_distances(self.prot2.positions[-1], verb_idxes)

        for value in [
            numverbs['VB'],  numverbs['VBD'],
            numverbs['VBG'], numverbs['VBN'],
            numverbs['VBP'], numverbs['VBZ'],
            maxscore, totalscore,
            int(cl1), int(far1),
            int(cl2), int(far2)]:
            if value >= 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(value)
            self.feat_current_col += 1

    def __keyword_count(self, mode="all"):

        keywords = dict({
        'Adenylate cyclase type 2': 0, 'Adenylate cyclase type 5': 0, 'Adenylate cyclase type 7': 0, 
        'Alpha-2C adrenergic receptor': 0, 'ATP-binding cassette sub-family F member 3': 0, 
        'Disintegrin and metalloproteinase domain-containing protein 2': 0,
        'Apoptotic chromatin condensation inducer in the nucleus': 0, 'Acetylcholine receptor subunit gamma': 0, 
        'Actin, cytoplasmic 1': 0, 'Actin-binding Rho-activating protein': 0, 'Beta-actin-like protein 2': 0, 
        'Actin-related protein T3': 0, 'Actin-like protein 6B': 0, 'Costars family protein ABRACL': 0
        })

        tokens = list()
        if mode == "all":
            tokens = self.prot1.sentence.tokens
        else:
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1
            tokens = self.prot1.sentence.tokens[init_coord:final_coord]

        for token in tokens:
            if token['lemma'] in keywords:
                keywords[token['lemma']] += 1

        for word, value in sorted(keywords.items()):
            if value > 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(value)
            self.feat_current_col += 1

    def features_todense(self):
        if self.features_sparse is None:
            self.compute_features()
        return self.features_sparse.todense()[0].tolist()[0]

    def predict(self):
        if self.features_sparse is None:
            self.compute_features()
        pred = InteractionCandidate.predictor.predict_proba(self.features_sparse)[:,1]

        self.votes = pred[0]
        self.__normalize_pred()
        if self.votes >= 0:
            self.label = True
        else:
            self.label = False


    def __normalize_pred(self):

        before = self.votes
        self.votes = (self.votes-0.55)/(1-0.55)
        self.votes = round(self.votes, 3)

        


    def to_html(self):

        init_coord  = self.between_idxes[0]
        final_coord = self.between_idxes[1]
        prot1_coords = [pos - 1 for pos in self.prot1.positions]
        prot2_coords = [pos - 1 for pos in self.prot2.positions]
        html_str = list()
        between = range(init_coord, final_coord)
        for i in range(0, len(self.prot1.sentence.tokens)):
            token = self.prot1.sentence.tokens[i]
            word = token['word']
            word = re.sub("-LRB-", "(", word)
            word = re.sub("-RRB-", ")", word)
            word = re.sub("-LSB-", "[", word)
            word = re.sub("-RSB-", "]", word)
            if i in between:
                if re.match('VB[DGNPZ]?', token['pos']):
                    # Verb in between
                    html_str.append('<span class="verb">%s</span>' % word)
                else:
                    # Normal word in between
                    html_str.append(word)
            elif i in prot1_coords:
                if i == prot1_coords[0]:
                    # Start of protein 1
                    html_str.append('<span class="prot"> %s' % word)
                    if i == prot1_coords[-1]:
                        # Protein of length 1
                        html_str.append('</span>')
                elif i == prot1_coords[-1]:
                    # End of protein of length > 1
                    html_str.append('%s </span>' % word)
                else:
                    # Middle of protein 1
                    html_str.append(word)
            elif i in prot2_coords:
                if i == prot2_coords[0]:
                    # Start of protein 2
                    html_str.append('<span class="prot"> %s' % word)
                    if i == prot2_coords[-1]:
                        # Protein of length 2
                        html_str.append('</span>')
                elif i == prot2_coords[-1]:
                    # End of protein of length > 2
                    html_str.append('%s </span>' % word)
                else:
                    # Middle of protein 2
                    html_str.append(word)
            else:
                # Neither verb nor protein
                html_str.append(word)
        return " ".join(html_str)



    def __str__(self):
        return "[%s] may interact with [%s]" % (self.prot1.symbol, self.prot2.symbol)


class TextNotAvailable(Exception):
    pass

class PubMedQueryError(Exception):
    pass

class ConnectionError(Exception):
    pass

class ProteinNotFound(Exception):
    pass

class GeneDictError(Exception):
    pass

class CannotLoadClassifier(Exception):
    pass
