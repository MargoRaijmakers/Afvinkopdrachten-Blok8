# Margo Raijmakers
# 12-04-2021
# Afvinkopdracht 2: Text Mining

import sys
from Bio import Entrez
from Bio import Medline
import matplotlib.pyplot as plt
import numpy as np


def ask_terms():
    """Deze functie vraagt op de zoektermen die gezocht moeten
    worden op PubMed.

    :return: de zoektermen in een lijst
    """
    terms = []
    try:
        terms_count = int(input("Aantal zoektermen: "))
        for i in range(terms_count):
            term = input("Zoekterm: ")
            terms.append(term)
    except ValueError:
        print("Niet een geldig getal.")
        sys.exit()
    return terms


def get_ids(term):
    """Deze functie haalt de ID's op van de PubMed artikelen waar de
    zoektermen in gevonden zijn.

    :param term: een term uit de lijst met opgegeven zoektermen
    :return: de ID's van de gevonden PubMed artikelen in een lijst
    """
    Entrez.email = "margoraijmakers@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=term,
                            retmax=100)  # maximaal 100 artikelen zodat
    # de grafiek sneller gemaakt is (kan aangepast worden)
    record = Entrez.read(handle)
    handle.close()
    idlist = record["IdList"]
    return idlist


def get_years(idlist):
    """Deze functie haalt de publicatiejaren op van de gevonden PubMed
    artikelen.

    :param idlist: de ID's van de gevonden PubMed artikelen in een lijst
    :return: de publicatiejaren in een lijst
    """
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",
                           retmode="text")
    records = Medline.parse(handle)
    records = list(records)
    years = []
    for record in records:
        try:
            years.append(int(record.get("DP").split(" ")[0]))
        except AttributeError:
            pass
        except ValueError:
            years.append(int(record.get("DP").split(" ")[1]))
    return years


def get_start_stop_year(years, pos):
    """Deze functie haalt de range van jaren op waarbinnen de PubMed
    artikelen gepubliceerd zijn.

    :param years: de publicatiejaren in een lijst
    :param pos: de positie in de lijst (0: startjaar, -1: stopjaar)
    :return: het jaar
    """
    if years[pos] % 10 > 5:
        year = years[pos] - years[pos] % 10 + 5
    else:
        if pos == -1:
            year = years[pos] - years[pos] % 10 + 5
        elif pos == 0 and years[0] % 10 == 0:
            year = years[pos] - 5
        else:
            year = years[pos] - years[pos] % 10
    return year


def make_year_count_dict(start_year, stop_year):
    """Deze functie maakt de dictionary aan voor het tellen van de
    hoeveelheid jaren per vijf jaar.

    :param start_year: het eerste jaar in de lijst met jaren afgerond
    op 5 of 0
    :param stop_year: het laatste jaar in de lijst met jaren afgerond
    op 5 of 0
    :return: de dictionary met de jaren opgesplitst op 5 jaar als keys
    met 0 als values
    """
    year_count = {}
    for i in range(start_year + 1, stop_year, 5):
        year_count["{}-{}".format(i, i + 4)] = 0
    return year_count


def count_years(years, year_count):
    """Deze functie telt het aantal jaren dat bij iedere vijf jaar
    hoort.

    :param years: de publicatiejaren in een lijst
    :param year_count: de dictionary met de jaren opgesplitst op 5 jaar
    als keys met 0 als values
    :return: de dictionary met de jaren opgesplitst op 5 jaar als keys
    met het aantal jaren als values
    """
    for year in years:
        for key, value in year_count.items():
            if int(key.split("-")[1]) >= year >= int(key.split("-")[0]):
                year_count[key] += 1
    return year_count


def get_unique_keys(year_counts):
    """Deze functie haalt de unieke keys op uit de year_counts
    dictionaries.

    :param year_counts: de dictionaries met de jaren opgesplitst op 5
    jaar als keys in een lijst
    met het aantal jaren als values
    :return: de unieke keys in een set
    """
    year_keys = []
    for year_count in year_counts:
        for key in year_count:
            year_keys.append(key)
    year_key_set = sorted(set(year_keys))
    return year_key_set


def hits_in_list(year_counts, year_key_set):
    """Deze functie plaatst de gevonden hits in een list.

    :param year_counts: de dictionaries met de jaren opgesplitst op 5
    jaar als keys in een lijst
    :param year_key_set: de unieke keys in een set
    :return: de hits in een lijst
    """
    hits_counts = []
    for year_count in year_counts:
        hits_count = []
        for year_key in year_key_set:
            # Als de key niet in de dictionary zit
            if year_count.get(year_key) is None:
                hits_count.append(0)
            # Als de key wel in de dictionary zit
            else:
                hits_count.append(year_count.get(year_key))
        hits_counts.append(hits_count[:])
        hits_count.clear()
    return hits_counts


def bar_plot(hits_counts, year_key_set, terms_with_results):
    """Deze functie maakt de grafiek van het aantal artikelen per 5
    jaar.

    :param hits_counts: de hits in een lijst
    :param year_key_set: de unieke keys in een set
    :param terms_with_results: de termen met resultaten
    """
    x_pos = [np.arange(len(hits_counts[0]))]
    for i in range(len(year_key_set)):
        x_pos.append([x + 0.25 for x in x_pos[i]])

    for i in range(len(hits_counts)):
        plt.bar(x_pos[i], hits_counts[i], width=0.25,
                label=terms_with_results[i])

    plt.ylabel('Hoeveelheid artikelen', fontweight='bold')
    plt.title('Hoeveelheid artikelen per vijf jaar', fontweight='bold')
    plt.xticks([r + 0.25 for r in range(len(hits_counts[0]))],
               year_key_set)
    plt.xticks(rotation=90, size=5.5)
    plt.legend()
    plt.show()


def main():
    terms = ask_terms()
    terms_with_results = []
    year_counts = []

    for term in terms:
        idlist = get_ids(term)
        years = get_years(idlist)
        # Als er geen resultaten gevonden zijn voor een zoekterm
        if len(years) == 0:
            print("Geen resultaten gevonden voor: {}".format(term))
        # Als er wel resultaten gevonden zijn voor een zoekterm
        else:
            terms_with_results.append(term)
            years.sort()
            start_year = get_start_stop_year(years, 0)
            stop_year = get_start_stop_year(years, -1)
            year_count = make_year_count_dict(start_year, stop_year)
            year_count = count_years(years, year_count)
            year_counts.append(year_count)

    # Als er geen resultaten zijn voor alle zoektermen
    if len(terms_with_results) == 0:
        sys.exit()
    # Als minimaal 1 zoekterm resultaten oplevert
    else:
        year_key_set = get_unique_keys(year_counts)
        hits_counts = hits_in_list(year_counts, year_key_set)
        bar_plot(hits_counts, year_key_set, terms_with_results)


main()
