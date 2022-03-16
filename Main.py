import numpy as np
import pandas
import os
import matplotlib.pyplot as plt
from Utility import Loader, Util, SignalProcessing, NoOffsetException


def find_correct_bound(lead_ecg, lead_ann):
    result = []
    for item in lead_ann:
        if type(item) is list:
            volt1 = lead_ecg[item[0]]
            volt2 = lead_ecg[item[1]]
            if abs(volt1) < abs(volt2):
                result.append(item[0])
            else:
                result.append(item[1])
        else:
            result.append(item)
    return result


def vote_among_leads(annotations, consensus_threshold: int):
    total = []
    for i in range(len(annotations)):
        ann_with_score = []
        for j in range(len(annotations[i])):
            score = 0
            ann_item = annotations[i][j]
            for k in range(len(annotations)):
                if k == i:
                    continue
                lead_ann = annotations[k]
                min_diff = Util.get_min_distance(ann_item, lead_ann)
                if min_diff < consensus_threshold:
                    score += 1
            ann_with_score.append((ann_item, score))
        total.append(ann_with_score)

    for i in range(len(total)):
        for j in range(len(total[i])):
            ann_item = total[i][j]
            score = ann_item[1]
            if score >= 8:
                continue
            else:
                high_conf_matches = []
                for k in range(len(total)):
                    if k == i:
                        continue
                    lead_ann = total[k]
                    match_ann = Util.get_closest_ann(ann_source=ann_item[0], input_list=lead_ann)
                    if match_ann[1] >= 8 and match_ann[1] != 12:
                        high_conf_matches.append(match_ann)
                if len(high_conf_matches) == 0:
                    # total[i][j] = (ann_item[0], -2)
                    v = 9
                else:
                    consensus_ann = round(
                        sum(ann_match for ann_match, score_match in high_conf_matches) / len(high_conf_matches))
                    if j - 1 >= 0 and abs(total[i][j - 1][0] - consensus_ann) <= consensus_threshold:
                        total[i][j] = (-1, -1)
                    elif j + 1 < len(total[i]) and abs(total[i][j + 1][0] - consensus_ann) <= consensus_threshold:
                        total[i][j] = (-1, -1)
                    else:
                        total[i][j] = (consensus_ann, 12)
    # Remove (-1, -1)
    # Union across 12 leads
    total_new = []
    for i in range(len(total)):
        row_ann = []
        for j in range(len(total[i])):
            ann_item = total[i][j]
            score = ann_item[1]
            if score == -1:
                continue
            else:
                row_ann.append(ann_item)
        total_new.append(row_ann)

    for i in range(len(total_new)):
        total_new[i].sort(key=lambda item: item[0])

    total_new_no_score = []
    for lead_ann in total_new:
        row_temp = []
        for item in lead_ann:
            row_temp.append(item[0])
        total_new_no_score.append(row_temp)

    return total_new, total_new_no_score


if __name__ == "__main__":
    ecg_ids = Util.get_ecg_list()
    column_names = ['P Start', 'P End', 'QRS Start', 'QRS End', 'T Start', 'T End']

    for ecg_id in ecg_ids:
        print('Processing ECG ' + str(ecg_id))
        ecg, frequency = Loader.get_ecg(ecg_id=ecg_id)
        ann = Loader.get_annotations(ecg_id=ecg_id)

        for index, row in ann.iterrows():
            lead_ecg = ecg[Util.get_lead_name(index=index)].values
            for col_name in column_names:
                lead_bound_candidates = row[col_name]
                corrected_bound = find_correct_bound(lead_ecg, lead_bound_candidates)
                ann.at[index, col_name] = corrected_bound

        consensus_threshold_qrs = round(frequency * 0.0625)
        consensus_threshold_t_p = round(consensus_threshold_qrs * 1.5)

        qrs_ons = ann['QRS Start'].values
        qrs_ends = ann['QRS End'].values
        qrs_ons_voted, qrs_ons = vote_among_leads(annotations=qrs_ons, consensus_threshold=consensus_threshold_qrs)
        qrs_ends_voted, qrs_ends = vote_among_leads(annotations=qrs_ends, consensus_threshold=consensus_threshold_qrs)

        t_ons = ann['T Start'].values
        t_ends = ann['T End'].values
        t_ons_voted, t_ons = vote_among_leads(annotations=t_ons, consensus_threshold=consensus_threshold_t_p)
        t_ends_voted, t_ends = vote_among_leads(annotations=t_ends, consensus_threshold=consensus_threshold_t_p)

        p_ons = ann['P Start'].values
        p_ends = ann['P End'].values
        p_ons_voted, p_ons = vote_among_leads(annotations=p_ons, consensus_threshold=consensus_threshold_t_p)
        p_ends_voted, p_ends = vote_among_leads(annotations=p_ends, consensus_threshold=consensus_threshold_t_p)

        # Preparing final annotations for saving
        for lead in range(12):
            ann.at[lead, 'P Start'] = p_ons[lead]
            ann.at[lead, 'P End'] = p_ends[lead]

            ann.at[lead, 'QRS Start'] = qrs_ons[lead]
            ann.at[lead, 'QRS End'] = qrs_ends[lead]

            ann.at[lead, 'T Start'] = t_ons[lead]
            ann.at[lead, 'T End'] = t_ends[lead]

        ann.insert(loc=0, column='Lead', value=['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'])
        ann.to_csv('Data/FinalAnnotation/' + str(ecg_id) + '_ann.csv', index=False)

        # Plotting QRS onset and offset and T-offset
        fig = plt.figure()
        fig.set_size_inches(18.5, 10.5)
        for lead in range(12):
            one_lead_signal = ecg.iloc[:, lead].values
            # one_lead_signal = one_lead_signal[:round(len(one_lead_signal) / 2) + 400]
            lead_onsets = qrs_ons[lead]
            lead_offsets = qrs_ends[lead]
            # lead_onsets.extend(t_ons[lead])
            lead_offsets.extend(t_ends[lead])
            my_plot = plt.subplot(6, 2, lead + 1)
            my_plot.scatter(lead_onsets, one_lead_signal[lead_onsets], c='r', marker="o", alpha=.5)
            my_plot.scatter(lead_offsets, one_lead_signal[lead_offsets], c='m', marker="o", alpha=.5)
            my_plot.plot(one_lead_signal)
        plt.show()
        # fig.savefig(str(ecg_id) + '.png', dpi=300)







